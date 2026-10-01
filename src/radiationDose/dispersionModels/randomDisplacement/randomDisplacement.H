/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::randomDisplacement

Description
    Random displacement model: turbulent dispersion as a diffusion of
    diffusivity K, with no velocity memory. Over a step of duration dt
    the particle is displaced by

        dx = (U + grad(K)) dt + sqrt(2 K dt) xi,

    xi a unit Gaussian triple drawn afresh every step: the Ito form of
    the diffusion equation, well mixed by construction (a uniform
    distribution of particles stays uniform however K varies). It is the
    diffusion limit of the well-mixed Langevin model (Wilson & Sawford
    1996), exact for times long against the eddy
    lifetime; it does not represent the ballistic spreading of a
    particle within one eddy.

    K is the discrete random walk's: K = sigma^2 tau_e / 2 = k tau_e / 3,
    with sigma = sqrt(2k/3) and tau_e = Cl k / epsilon or Cl / (Cmu omega)
    (see eddyDiffusivity), so the two models spread particles at the same
    rate over long times. k and the dissipation are interpolated linearly
    in the particle's tetrahedron and grad(K) is the exact gradient of the
    K so formed.

    Where the discrete random walk holds a velocity for an eddy lifetime,
    an eddy can be as long as a third of the distance to a wall in a
    resolved boundary layer, and a held eddy carries particles into the
    viscous sublayer faster than the local turbulence would; the drift
    correction is only first order in that ratio. This model has no
    eddies to hold.

    The CFL bound sees only the deterministic part U + grad(K): bounding
    the step by the random draw it is about to take would make large
    draws take short steps and under-spread the particles. Instead the
    step is limited so that the diffusive displacement sqrt(2 K dt) does
    not exceed the CFL displacement (cflMax times the cell size).

    References:
        Wilson, J. D., & Sawford, B. L. (1996). Review of Lagrangian
        stochastic models for trajectories in the turbulent atmosphere.
        Boundary-Layer Meteorology, 78, 191-210.

Usage
    dispersion
    {
        type        randomDisplacement;
        k           k;            // optional, default "k"
        omega       omega;        // or epsilon (default "epsilon")
        Cmu         0.09;         // optional, default 0.09 (with omega)
        Cl          0.15;         // optional, default 0.15
    }

\*---------------------------------------------------------------------------*/

#ifndef dose_randomDisplacement_H
#define dose_randomDisplacement_H

#include "dispersionModel.H"
#include "eddyDiffusivity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                      Class randomDisplacement Declaration
\*---------------------------------------------------------------------------*/

class randomDisplacement
:
    public dispersionModel
{
public:

    // Public Classes

        //- Per-track state: the current step's sample and draw
        class RDState : public dispersionModel::State
        {
        public:

            //- Diffusivity at the step's start [m^2/s]
            scalar K;

            //- Drift velocity grad(K) at the step's start [m/s]
            vector drift;

            //- Unit Gaussian draw for the step
            vector xi;

            RDState()
            :
                K(0),
                drift(vector::zero),
                xi(vector::zero)
            {}
        };


private:

    // Private Data

        //- k, the dissipation, tau_e and K = k tau_e / 3
        eddyDiffusivity turbulence_;


public:

    //- Runtime type information
    TypeName("randomDisplacement");


    // Constructors

        randomDisplacement(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~randomDisplacement() = default;


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
            return autoPtr<dispersionModel::State>(new RDState());
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

        //- grad(K) + sqrt(2 K / dt) xi: the step's displacement over dt
        virtual vector fluctuation
        (
            const dispersionModel::State& state,
            scalar dt
        ) const;

        //- grad(K) alone: the CFL bound must not see the draw
        virtual vector cflFluctuation
        (
            const dispersionModel::State& state,
            scalar dt
        ) const;

        //- maxDisplacement^2 / (2 K): sqrt(2 K dt) stays within it
        virtual scalar maxStep
        (
            const dispersionModel::State& state,
            scalar maxDisplacement
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
