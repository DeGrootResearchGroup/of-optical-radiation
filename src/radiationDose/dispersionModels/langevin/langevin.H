/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::langevin

Description
    Well-mixed Langevin model of the fluctuating velocity (Thomson 1987),
    for isotropic Gaussian turbulence whose variance sigma^2 = 2k/3 varies
    in space. The particle carries a normalized velocity v = u'/sigma that
    relaxes on the Lagrangian time scale T_L:

        dv = (-v / T_L + grad(sigma)) dt + sqrt(2 / T_L) dW,
        u' = sigma v.

    This is the unique first-order model that keeps a uniform distribution
    of particles uniform in such turbulence: the drift grad(sigma) on v
    is the well-mixed condition's, and u' changes continuously as the
    particle moves into weaker or stronger turbulence. It keeps the
    velocity memory the discrete random walk has (u' is correlated over
    about T_L, so a particle spreads as sigma^2 t^2 at first and as
    2 K t once t >> T_L) without holding a velocity through an eddy whose
    length is comparable with the distance over which the turbulence
    changes, which near a resolved wall carries particles into the viscous
    sublayer.

    T_L = tau_e / 2, with tau_e the discrete random walk's eddy lifetime
    (Cl k / epsilon or Cl / (Cmu omega), see eddyDiffusivity), so the
    long-time diffusivity sigma^2 T_L is the walk's K = k tau_e / 3 and
    the three dispersion models spread particles at the same rate over
    long times.

    Each step integrates v and the displacement exactly with sigma,
    grad(sigma) and T_L frozen at the step's start -- the end velocity and
    the displacement are jointly Gaussian with known means, variances and
    covariance -- and adds the mean displacement that freezing them
    leaves out, from expanding sigma and T_L to first order along the
    path with v at its stationary statistics:

        dt [psi1(h) sigma T_L grad(sigma) + psi2(h) sigma^2 grad(T_L)],
        psi1 = 1 - (1 - e^-h)/h,   psi2 = 1 - 2 (1 - e^-h)/h + e^-h,

    with h = dt / T_L. Both vanish for h << 1 (psi1 ~ h/2, psi2 ~ h^2/6),
    where the step is the Langevin model's, and tend to one for h >> 1,
    where the mean drift becomes sigma T_L grad(sigma) + sigma T_L
    grad(sigma) + sigma^2 grad(T_L) = grad(K) and the step is the random
    displacement model's with the same K. One integrator therefore covers
    both limits and every step between them, with no switch between
    models.

    The step is at most maxStepFraction max(T_L, minLagrangianTime), and
    at most maxDisplacement^2 / (2 K) so that the random spread of a step
    stays within the distance over which the first-order expansion holds.
    It is continuous in space: a step that jumps between neighbouring
    positions unbalances the particle fluxes across the jump, and a step
    that varies with position along a path whose velocity is still
    remembered drifts particles towards shorter steps by about
    -sigma^2 grad(dt) / 2 (Wilson and Flesch 1993), which a small
    maxStepFraction keeps small. Where T_L is below minLagrangianTime --
    next to a resolved wall, where k-omega SST gives microseconds -- the
    step stops shrinking, and h grows towards the diffusion limit.

    The CFL bound sees only the deterministic part of the step's motion,
    never its random draws (a bound by the draw it is about to take would
    give large draws short steps and bias the spreading). At walls the
    velocity is reflected with the particle.

    References:
        Thomson, D. J. (1987). Criteria for the selection of stochastic
        models of particle trajectories in turbulent flows. Journal of
        Fluid Mechanics, 180, 529-556.

        Wilson, J. D., & Flesch, T. K. (1993). Flow boundaries in
        random-flight dispersion models: enforcing the well-mixed
        condition. Journal of Applied Meteorology, 32, 1695-1707.

        Wilson, J. D., & Sawford, B. L. (1996). Review of Lagrangian
        stochastic models for trajectories in the turbulent atmosphere.
        Boundary-Layer Meteorology, 78, 191-210.

Usage
    dispersion
    {
        type                langevin;
        k                   k;        // optional, default "k"
        omega               omega;    // or epsilon (default "epsilon")
        Cmu                 0.09;     // optional, default 0.09 (with omega)
        Cl                  0.15;     // optional, default 0.15
        maxStepFraction     0.02;     // optional: step <= this x T_L
        minLagrangianTime   1e-2;     // optional [s]: the step's T_L floor
    }

\*---------------------------------------------------------------------------*/

#ifndef dose_langevin_H
#define dose_langevin_H

#include "dispersionModel.H"
#include "eddyDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                         Class langevin Declaration
\*---------------------------------------------------------------------------*/

class langevin
:
    public dispersionModel
{
public:

    // Public Classes

        //- Per-track state
        class LangevinState : public dispersionModel::State
        {
        public:

            //- Normalized velocity u'/sigma
            vector v;

            //- False until v has been drawn from its stationary
            //  distribution at the track's first step
            bool hasVelocity;

            // The current step, set by beginStep()

                //- sigma [m/s], grad(sigma) [1/s], T_L [s] and grad(T_L)
                //  [s/m] at the step's start
                scalar sigma;
                vector gradSigma;
                scalar TL;
                vector gradTL;

                //- K [m^2/s] at the step's start
                scalar K;

                //- Unit Gaussian draws for the step
                vector xi1;
                vector xi2;

            LangevinState()
            :
                v(vector::zero),
                hasVelocity(false),
                sigma(0),
                gradSigma(vector::zero),
                TL(0),
                gradTL(vector::zero),
                K(0),
                xi1(vector::zero),
                xi2(vector::zero)
            {}
        };


private:

    // Private Data

        //- k, the dissipation, tau_e and K = k tau_e / 3
        eddyDiffusivity turbulence_;

        //- Largest step as a fraction of T_L
        scalar maxStepFraction_;

        //- Below this T_L [s] the step stops shrinking
        scalar minLagrangianTime_;


    // Private Member Functions

        //- The step's update over dt: the mean of the normalized
        //  velocity's time integral (Xm, including the first-order
        //  correction for sigma and T_L varying along the path), its
        //  random part (Xn) and the end velocity (vEnd)
        void integrate
        (
            const LangevinState& s,
            scalar dt,
            vector& Xm,
            vector& Xn,
            vector& vEnd
        ) const;


public:

    //- Runtime type information
    TypeName("langevin");


    // Constructors

        langevin(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~langevin() = default;


    // Member Functions

        //- Static requiredFields: k and the dissipation the dict names
        static wordList requiredFields(const dictionary& dict)
        {
            return eddyDiffusivity::requiredFields(dict);
        }

        virtual wordList requiredFields() const
        {
            return turbulence_.requiredFields();
        }

        virtual autoPtr<dispersionModel::State> newState() const
        {
            return autoPtr<dispersionModel::State>(new LangevinState());
        }

        //- Build the k and dissipation interpolators
        virtual void correct();

        virtual void beginStep
        (
            dispersionModel::State& state,
            const barycentric& coordinates,
            const tetIndices& tetIs,
            scalar dtMax,
            randomGenerator& rng
        ) const;

        //- The mean u' over the step: its displacement over dt
        virtual vector fluctuation
        (
            const dispersionModel::State& state,
            scalar dt
        ) const;

        //- The deterministic part of fluctuation()
        virtual vector cflFluctuation
        (
            const dispersionModel::State& state,
            scalar dt
        ) const;

        //- The smaller of maxStepFraction max(T_L, minLagrangianTime)
        //  and maxDisplacement^2 / (2 K)
        virtual scalar maxStep
        (
            const dispersionModel::State& state,
            scalar maxDisplacement
        ) const;

        //- Advance v by the step
        virtual void endStep(dispersionModel::State& state, scalar dt) const;

        //- Mirror v
        virtual void reflect
        (
            dispersionModel::State& state,
            const vector& n
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
