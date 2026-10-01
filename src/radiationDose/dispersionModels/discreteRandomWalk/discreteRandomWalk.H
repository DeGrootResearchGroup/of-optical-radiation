/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::discreteRandomWalk

Description
    Discrete Random Walk (DRW) dispersion model in the Gosman-Ioannides
    sense. The fluctuation velocity u' is held constant through an eddy
    lifetime tau_e = C_L * k / epsilon, then resampled. Each component
    of u' is drawn from N(0, sigma) with sigma^2 = 2k/3 (isotropic-
    turbulence assumption).

    The dissipation comes from the turbulence model's own field: the
    rate epsilon of a k-epsilon model, or the specific rate omega of a
    k-omega model, with epsilon = Cmu * k * omega. Given omega, the
    lifetime is tau_e = C_L / (Cmu * omega), the same eddy as
    C_L * k / epsilon. Name one of the two; epsilon is the default.
    The fields are read through eddyDiffusivity, shared with the
    randomDisplacement model.

    A walk of eddies held for tau_e spreads particles with the
    diffusivity K = sigma^2 tau_e / 2 = k tau_e / 3. Where K varies in
    space, a walk that only adds u' violates the well-mixed condition
    (Thomson 1987): a uniform distribution of particles does not stay
    uniform, and particles collect where K is small. In the viscous
    sublayer of a mesh resolved to the wall, where k and the mean
    velocity both vanish, they enter the first cells and barely leave.

    With wellMixed on (the default) the walk is corrected in three
    ways, which together keep a uniform distribution uniform:

    1. Drift. The drift velocity grad(K) is added to u'. It is the
       Ito drift that makes a walk of diffusivity K well mixed: the
       inhomogeneous-turbulence drift of Legg & Raupach (1982) and
       MacInnes & Bracco (1992), tau grad(sigma^2), with the gradient
       of tau_e kept as well. k and the dissipation are interpolated
       linearly in the particle's tetrahedron, and the drift is the
       exact gradient of the K so interpolated -- the K the eddies
       actually realise -- so the two cannot disagree across cells.

    2. Exact eddy accounting. The displacement over an outer step of
       duration dt is the exact integral of the piecewise-constant
       eddy velocity: the rest of the current eddy, then every whole
       eddy that fits (the sum of m of them is one Gaussian of
       variance m sigma^2 tau_e^2), then the start of the eddy that
       carries into the next step. The realised diffusivity is then
       k tau_e / 3 whatever dt is. Without it, wherever tau_e < dt the
       walk resamples once per step and holds u' for the whole step,
       so its diffusivity is sigma^2 dt / 2: set by the step, not by
       the turbulence. On a wall-resolved k-omega SST mesh tau_e is
       microseconds near the wall, and shorter than a millisecond-scale
       dtMax in most of the volume.

    3. Reflection. When the particle reflects off a wall, the eddy
       velocity it carries is reflected too. Otherwise an eddy aimed
       at the wall keeps pressing the particle against it until the
       eddy expires.

    With wellMixed off the walk is the uncorrected one: u' sampled
    from the owning cell's k and dissipation, resampled once the eddy
    has used up its lifetime in steps of dtMax, no drift and no
    reflection of u'.

    Per-particle eddy memory (the current u', its remaining lifetime
    and the step's draws) lives on the track object as a DRWState. The
    model has no shared mutable state across tracks, so the per-track
    integration loop is thread-parallel; the interpolators built by
    correct() are only read during it.

    References:
        Thomson, D. J. (1987). Criteria for the selection of
        stochastic models of particle trajectories in turbulent flows.
        Journal of Fluid Mechanics, 180, 529-556.

        Legg, B. J., & Raupach, M. R. (1982). Markov-chain simulation
        of particle dispersion in inhomogeneous flows: the mean drift
        velocity induced by a gradient in Eulerian velocity variance.
        Boundary-Layer Meteorology, 24, 3-13.

        MacInnes, J. M., & Bracco, F. V. (1992). Stochastic particle
        dispersion modeling and the tracer-particle limit. Physics of
        Fluids A, 4(12), 2809-2824.

Usage
    dispersion
    {
        type        discreteRandomWalk;
        k           k;            // optional, default "k"
        epsilon     epsilon;      // optional, default "epsilon"
        Cl          0.15;         // optional, default 0.15
        wellMixed   true;         // optional, default true
    }

    or, for a k-omega model (kOmegaSST, say):

    dispersion
    {
        type        discreteRandomWalk;
        omega       omega;        // instead of epsilon; not both
        Cmu         0.09;         // optional, default 0.09
    }

\*---------------------------------------------------------------------------*/

#ifndef dose_discreteRandomWalk_H
#define dose_discreteRandomWalk_H

#include "dispersionModel.H"
#include "eddyDiffusivity.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class discreteRandomWalk Declaration
\*---------------------------------------------------------------------------*/

class discreteRandomWalk
:
    public dispersionModel
{
public:

    // Public Classes

        //- Per-track eddy state. Lives on the track object.
        class DRWState : public dispersionModel::State
        {
        public:

            //- Velocity of the current eddy [m/s]
            vector uPrime;

            //- Time the current eddy has left [s]
            scalar remaining;

            // The current step (well-mixed walk), set by beginStep()

                //- sqrt(2k/3) [m/s] and tau_e [s] at the step's start,
                //  for the eddies that begin during the step
                scalar sigma;
                scalar tau;

                //- Drift velocity grad(K) at the step's start [m/s]
                vector drift;

                //- Unit Gaussian draws: the sum of the step's whole
                //  eddies, and the eddy that carries over
                vector xiSum;
                vector xiNew;

                //- The draws are waiting to be used. A step the CFL
                //  bound shortens may end before the current eddy does
                //  and leave them unused; they are then kept for the
                //  next step, not redrawn. Redrawing would keep only
                //  the draws that gave a short enough step -- those
                //  pointing against the current eddy -- and bias the
                //  next eddy towards reversing it.
                bool drawn;

            DRWState()
            :
                uPrime(vector::zero),
                remaining(0),
                sigma(0),
                tau(0),
                drift(vector::zero),
                xiSum(vector::zero),
                xiNew(vector::zero),
                drawn(false)
            {}
        };


private:

    // Private Data

        //- k, the dissipation, tau_e and K = k tau_e / 3
        eddyDiffusivity turbulence_;

        //- Apply the drift, exact eddy accounting and reflection that
        //  keep the walk well mixed (see Description)
        Switch wellMixed_;


    // Private Member Functions

        //- Uncorrected walk: resample from the owning cell's fields
        //  once the eddy has expired, then age it by dtMax
        void beginStepUncorrected
        (
            DRWState& s,
            label celli,
            scalar dtMax,
            randomGenerator& rng
        ) const;

        //- For a step of duration dt that outlasts the current eddy:
        //  the number m of whole new eddies that fit in it, and the
        //  time `frac` the eddy that carries over has already run.
        //  False when the step's turbulence is nil (no new eddies).
        bool newEddies
        (
            const DRWState& s,
            scalar dt,
            scalar& m,
            scalar& frac
        ) const;


public:

    //- Runtime type information
    TypeName("discreteRandomWalk");


    // Constructors

        discreteRandomWalk(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~discreteRandomWalk() = default;


    // Member Functions

        //- Static requiredFields: k and the dissipation the dict
        //  names. Registered in the base class's dispatch table at
        //  static-init time.
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
            return autoPtr<dispersionModel::State>(new DRWState());
        }

        //- Build the k and dissipation interpolators from the
        //  registry's current fields (well-mixed walk)
        virtual void correct();

        virtual void beginStep
        (
            dispersionModel::State& state,
            const barycentric& coordinates,
            const tetIndices& tetIs,
            scalar dtMax,
            randomGenerator& rng
        ) const;

        virtual vector fluctuation
        (
            const dispersionModel::State& state,
            scalar dt
        ) const;

        virtual void endStep(dispersionModel::State& state, scalar dt) const;

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
