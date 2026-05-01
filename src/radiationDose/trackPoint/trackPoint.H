/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

Class
    Foam::dose::trackPoint

Description
    One vertex of a particle trajectory.

    Position is in mesh-native units (m, OpenFOAM convention). Time is in
    seconds. Accumulated dose is in mJ/cm^2 (UV convention; see radiationDose
    function object for the unit-conversion convention).

\*---------------------------------------------------------------------------*/

#ifndef trackPoint_H
#define trackPoint_H

#include "vector.H"
#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                         Struct trackPoint Declaration
\*---------------------------------------------------------------------------*/

struct trackPoint
{
    //- Position [m]
    vector x;

    //- Time since seeding [s]
    scalar t;

    //- Accumulated dose [mJ/cm^2]
    scalar D;

    //- Last known cell index (-1 if unset)
    label celli;

    trackPoint()
    :
        x(vector::zero),
        t(0),
        D(0),
        celli(-1)
    {}

    trackPoint(const vector& x_, scalar t_, scalar D_, label celli_)
    :
        x(x_),
        t(t_),
        D(D_),
        celli(celli_)
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
