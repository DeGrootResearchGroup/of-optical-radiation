/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "schillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(schillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, schillerNaumann, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::schillerNaumann::schillerNaumann(const dictionary& /*dict*/)
:
    dragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dose::schillerNaumann::tauP
(
    const vector& V,
    const vector& Useen,
    scalar rhoP,
    scalar dP,
    scalar nuF,
    scalar muF
) const
{
    const scalar tauStokes = rhoP*dP*dP/(18.0*muF);
    const scalar Rep = mag(Useen - V)*dP/max(nuF, small);
    const scalar correction = 1.0 + 0.15*pow(Rep, 0.687);
    return tauStokes/correction;
}


// ************************************************************************* //
