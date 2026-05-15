/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::trackPoint

Description
    One vertex of a particle trajectory.

    Position is stored single-precision (Vector<float>) -- the only
    consumer is the VTK writer, which emits POINTS as float anyway, and
    the per-particle DynamicList<trackPoint> is the dominant memory
    term for long-residence runs. Time and dose stay double-precision
    because they are integrated into double-precision accumulators on
    the particle; the trackPoint just holds the per-vertex snapshot.

    Position is in mesh-native units (m, OpenFOAM convention). Time is
    in seconds. Accumulated dose is in mJ/cm^2 (UV convention; see
    radiationDose function object for the unit-conversion convention).

\*---------------------------------------------------------------------------*/

#ifndef dose_trackPoint_H
#define dose_trackPoint_H

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
    //- Position [m], single-precision
    Vector<float> x;

    //- Time since seeding [s]
    scalar t;

    //- Accumulated dose [mJ/cm^2]
    scalar D;

    //- Last known cell index (-1 if unset)
    label celli;

    trackPoint()
    :
        x(0.0f, 0.0f, 0.0f),
        t(0),
        D(0),
        celli(-1)
    {}

    trackPoint(const vector& x_, scalar t_, scalar D_, label celli_)
    :
        x(float(x_.x()), float(x_.y()), float(x_.z())),
        t(t_),
        D(D_),
        celli(celli_)
    {}

    //- Cartesian-vector view of the stored position. Used by the
    //  seeding -> meshSearch path that needs a double-precision
    //  position to locate the surrounding tet, and by the VTK gather
    //  path that flattens vertices into a List<vector>.
    inline vector xd() const
    {
        return vector(x.x(), x.y(), x.z());
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
