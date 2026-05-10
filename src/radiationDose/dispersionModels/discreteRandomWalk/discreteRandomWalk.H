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
    sense. The fluctuation velocity u' is held constant through a
    randomly drawn eddy lifetime tau_e = C_L * k / epsilon, then
    resampled. Each component of u' is drawn from N(0, sigma) with
    sigma^2 = 2k/3 (isotropic-turbulence assumption).

    Per-particle eddy memory (current u' and remaining lifetime) lives
    on the track object as a DRWState. The model has no shared mutable
    state across tracks, which makes the per-track integration loop
    trivially thread-parallel.

Usage
    dispersion
    {
        type        discreteRandomWalk;
        k           k;            // optional, default "k"
        epsilon     epsilon;      // optional, default "epsilon"
        Cl          0.15;         // optional, default 0.15
    }

\*---------------------------------------------------------------------------*/

#ifndef dose_discreteRandomWalk_H
#define dose_discreteRandomWalk_H

#include "dispersionModel.H"

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
    // Private Data

        //- Name of the turbulent kinetic energy field [m^2/s^2]
        word kName_;

        //- Name of the turbulent dissipation rate field [m^2/s^3]
        word epsilonName_;

        //- Eddy lifetime constant: tau_e = Cl * k / epsilon
        scalar Cl_;


public:

    // Public Classes

        //- Per-track eddy state. Lives on the track object.
        class DRWState : public dispersionModel::State
        {
        public:
            vector uPrime;
            scalar remaining;

            DRWState()
            :
                uPrime(vector::zero),
                remaining(0)
            {}
        };


    //- Runtime type information
    TypeName("discreteRandomWalk");


    // Constructors

        discreteRandomWalk(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~discreteRandomWalk() = default;


    // Member Functions

        virtual wordList requiredFields() const
        {
            return wordList({kName_, epsilonName_});
        }

        virtual autoPtr<dispersionModel::State> newState() const
        {
            return autoPtr<dispersionModel::State>(new DRWState());
        }

        virtual vector fluctuation
        (
            dispersionModel::State& state,
            const vector& x,
            label celli,
            scalar dt,
            randomGenerator& rng
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
