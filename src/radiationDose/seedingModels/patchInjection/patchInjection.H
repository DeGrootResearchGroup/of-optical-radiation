/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::patchInjection

Description
    Seeds particles uniformly across one or more boundary patches,
    weighted by face area. Each particle's initial position is sampled
    randomly within a face (fan triangulation around the face centroid),
    and the initial cell index is taken as the face's owner cell.

Usage
    seeding
    {
        type        patchInjection;
        patches     (inlet);
        nParticles  10000;
    }

\*---------------------------------------------------------------------------*/

#ifndef dose_patchInjection_H
#define dose_patchInjection_H

#include "seedingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class patchInjection Declaration
\*---------------------------------------------------------------------------*/

class patchInjection
:
    public seedingModel
{
    // Private Data

        //- List of boundary patch names to seed from
        wordList patchNames_;

        //- Total number of particles to launch
        label nParticles_;


public:

    //- Runtime type information
    TypeName("patchInjection");


    // Constructors

        patchInjection(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~patchInjection() = default;


    // Member Functions

        virtual List<trackPoint> seed(randomGenerator& rng) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
