/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::dispersionModel

Description
    Abstract base class for stochastic turbulent dispersion models.
    Returns a velocity fluctuation u' to add to the mean velocity at
    a given particle position.

\*---------------------------------------------------------------------------*/

#ifndef dose_dispersionModel_H
#define dose_dispersionModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "randomGenerator.H"
#include "runTimeSelectionTables.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class dispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersionModel
{

protected:

    // Protected Data

        //- Model dictionary (the contents of the "dispersion" sub-dict)
        const dictionary dict_;

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("dispersionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dispersionModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        dispersionModel(const dictionary& dict, const fvMesh& mesh);


    //- Selector
    static autoPtr<dispersionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~dispersionModel();


    // Member Functions

        //- Names of fields the model needs loaded into the registry
        //  (e.g. k, epsilon for DRW). Default: none.
        virtual wordList requiredFields() const
        {
            return wordList();
        }

        //- Reset any per-track internal state (e.g. eddy memory).
        //  Called by the function object at the start of every execute().
        virtual void reset()
        {}

        //- Velocity fluctuation u' [m/s] at position x in cell celli over a
        //  step of size dt [s]. trackId disambiguates per-particle state for
        //  models that maintain memory across substeps (e.g. DRW eddies).
        //  Concrete models look up their own turbulence fields (k, epsilon,
        //  ...) from the registry as needed.
        virtual vector fluctuation
        (
            label trackId,
            const vector& x,
            label celli,
            scalar dt,
            randomGenerator& rng
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
