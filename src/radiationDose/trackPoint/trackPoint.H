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

    Position is in mesh-native units (m, OpenFOAM convention). Time is in
    seconds. Accumulated dose is in mJ/cm^2 (UV convention; see radiationDose
    function object for the unit-conversion convention).

\*---------------------------------------------------------------------------*/

#ifndef dose_trackPoint_H
#define dose_trackPoint_H

#include "vector.H"
#include "scalar.H"
#include "label.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"

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


// Required by Foam::UList<T>::operator>> for uniform-list detection
// during stream parsing. Field-wise comparison is the natural fit for
// this plain-data vertex record.
inline bool operator==(const trackPoint& a, const trackPoint& b)
{
    return a.x == b.x && a.t == b.t && a.D == b.D && a.celli == b.celli;
}

inline bool operator!=(const trackPoint& a, const trackPoint& b)
{
    return !(a == b);
}


// Stream operators. trackPoint is the element type of the
// DynamicList<trackPoint> trajectory carried on each dosePathParticle;
// when the particle is serialised across a processor patch (Pstream),
// the list invokes these element operators in turn so the full pre-
// handoff trajectory survives the cross-rank transfer.
inline Ostream& operator<<(Ostream& os, const trackPoint& p)
{
    os  << p.x << token::SPACE
        << p.t << token::SPACE
        << p.D << token::SPACE
        << p.celli;
    return os;
}

inline Istream& operator>>(Istream& is, trackPoint& p)
{
    is >> p.x >> p.t >> p.D >> p.celli;
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
