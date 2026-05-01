/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::seedingModel

Description
    Abstract base class for particle seeding strategies. A seedingModel
    produces an initial list of trackPoints (one per particle) when the
    radiationDose function object's execute() is invoked.

\*---------------------------------------------------------------------------*/

#ifndef dose_seedingModel_H
#define dose_seedingModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "randomGenerator.H"
#include "runTimeSelectionTables.H"
#include "trackPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class seedingModel Declaration
\*---------------------------------------------------------------------------*/

class seedingModel
{

protected:

    // Protected Data

        //- Model dictionary (the contents of the "seeding" sub-dict)
        const dictionary dict_;

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("seedingModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        seedingModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from dictionary and mesh
        seedingModel(const dictionary& dict, const fvMesh& mesh);


    //- Selector
    static autoPtr<seedingModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~seedingModel();


    // Member Functions

        //- Generate initial trackPoints (one per particle to launch).
        //  rng is the shared RNG so that the function-object's "seed"
        //  reproduces deterministically across stages.
        virtual List<trackPoint> seed(randomGenerator& rng) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
