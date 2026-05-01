/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "track.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

const Foam::NamedEnum<Foam::dose::track::endReason, 6>
Foam::dose::track::endReasonNames
{
    "active",
    "escaped",
    "timedOut",
    "stuck",
    "leftDomain",
    "terminated"
};


// ************************************************************************* //
