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

    Each step integrates v and the displacement exactly with the
    coefficients (sigma, grad(sigma), T_L) frozen at the step's start:
    the end velocity and the displacement are jointly Gaussian with known
    means, variances and covariance. Frozen coefficients are consistent
    only on steps short against T_L, so the step is limited to
    maxStepFraction T_L. Where T_L is shorter than minLagrangianTime --
    next to a resolved wall, where k-omega SST gives microseconds -- the
    step instead follows the diffusion limit with the same K,
    dx = grad(K) dt + sqrt(2 K dt) xi, and v is drawn afresh from its
    stationary distribution N(0, 1) when the particle leaves that zone.
    Eddies there are far shorter than any step, so the two descriptions
    agree where they meet.

    The CFL bound sees only the deterministic part of the step's motion,
    never its random draws (a bound by the draw it is about to take would
    give large draws short steps and bias the spreading). At walls the
    velocity is reflected with the particle.

    References:
        Thomson, D. J. (1987). Criteria for the selection of stochastic
        models of particle trajectories in turbulent flows. Journal of
        Fluid Mechanics, 180, 529-556.

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
        maxStepFraction     0.05;     // optional: step <= this x T_L
        minLagrangianTime   1e-4;     // optional [s]: diffusion below it
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

            //- False until v has been drawn: a new track, or one that
            //  has been in the diffusion zone
            bool hasVelocity;

            // The current step, set by beginStep()

                //- True when the step follows the diffusion limit
                bool diffusive;

                //- sigma [m/s], grad(sigma) [1/s] and T_L [s] at the
                //  step's start
                scalar sigma;
                vector gradSigma;
                scalar TL;

                //- K [m^2/s] and grad(K) [m/s] at the step's start
                scalar K;
                vector gradK;

                //- Unit Gaussian draws for the step
                vector xi1;
                vector xi2;

            LangevinState()
            :
                v(vector::zero),
                hasVelocity(false),
                diffusive(false),
                sigma(0),
                gradSigma(vector::zero),
                TL(0),
                K(0),
                gradK(vector::zero),
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

        //- Below this T_L [s] the step follows the diffusion limit
        scalar minLagrangianTime_;


    // Private Member Functions

        //- The step's exact Langevin update over dt: the mean of the
        //  normalized velocity's time integral (Xm), its random part
        //  (Xn) and the end velocity (vEnd)
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

        //- maxStepFraction T_L, or in the diffusion zone
        //  maxDisplacement^2 / (2 K)
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
