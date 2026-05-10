/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::pointInjection

Description
    Seeds particles uniformly inside a small region centred at a user-
    specified point in the interior of the domain. Region shape is
    either a sphere (`type sphere; radius R`) or an axis-aligned box
    (`type box; size (Lx Ly Lz)`, full edge lengths). Points are drawn
    by rejection sampling in the region's bounding cube and accepted
    only if they fall inside both the shape and a local mesh cell;
    candidates that miss the shape, or that map to no cell on the
    current rank, are redrawn.

    For the typical case (region fully interior, larger than a single
    cell) rejection sampling is exactly uniform in space and the only
    overhead is at most ~1.9 candidate draws per particle (sphere; box
    has 100 % acceptance). The maxAttempts guard prevents infinite
    loops when the region is entirely outside the local mesh.

Usage
    seeding
    {
        type        pointInjection;
        nParticles  10000;
        region
        {
            type    sphere;
            centre  (0.5 0 0);
            radius  0.05;
        }
        // or:
        // region
        // {
        //     type    box;
        //     centre  (0.5 0 0);
        //     size    (0.1 0.05 0.05);
        // }
    }

\*---------------------------------------------------------------------------*/

#ifndef dose_pointInjection_H
#define dose_pointInjection_H

#include "seedingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class pointInjection Declaration
\*---------------------------------------------------------------------------*/

class pointInjection
:
    public seedingModel
{
public:

    enum class shapeType { sphere, box };


private:

    // Private Data

        //- Total number of particles to launch (across all ranks)
        label nParticles_;

        //- Region shape
        shapeType shape_;

        //- Region centre
        point centre_;

        //- Sphere radius (used when shape_ == sphere)
        scalar radius_;

        //- Box full edge lengths (used when shape_ == box)
        vector size_;

        //- Per-particle attempt cap (safeguard against infinite loops
        //  when the region is mostly or entirely outside the local mesh)
        label maxAttempts_;


public:

    //- Runtime type information
    TypeName("pointInjection");


    // Constructors

        pointInjection(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~pointInjection() = default;


    // Member Functions

        virtual List<trackPoint> seed(randomGenerator& rng) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
