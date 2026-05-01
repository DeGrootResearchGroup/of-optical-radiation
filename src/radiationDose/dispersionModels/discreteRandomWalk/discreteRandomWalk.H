/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::discreteRandomWalk

Description
    Discrete Random Walk (DRW) dispersion model in the Gosman-Ioannides
    sense. The fluctuation velocity u' is held constant through a
    randomly drawn eddy lifetime tau_e = C_L * k / epsilon, then
    resampled. Each component of u' is drawn from N(0, sigma) with
    sigma^2 = 2k/3 (isotropic-turbulence assumption).

    Per-particle eddy memory (current u', remaining lifetime) is held
    in a Map keyed by track ID. The Map is cleared by reset() at the
    start of every execute() call so successive runs of the same
    function object don't carry stale state.

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

#include "Map.H"
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

        //- Per-track eddy state. Mutable so fluctuation() can stay const.
        struct eddyState
        {
            vector uPrime;
            scalar remaining;
        };
        mutable Map<eddyState> states_;


public:

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

        virtual void reset();

        virtual vector fluctuation
        (
            label trackId,
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
