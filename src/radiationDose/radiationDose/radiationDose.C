/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "radiationDose.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(radiationDose, 0);
    addToRunTimeSelectionTable(functionObject, radiationDose, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::radiationDose::radiationDose
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_("U"),
    GName_("G")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::functionObjects::radiationDose::~radiationDose()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::radiationDose::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    UName_ = dict.lookupOrDefault<word>("U", "U");
    GName_ = dict.lookupOrDefault<word>("fluenceRate", "G");

    return true;
}


Foam::wordList Foam::functionObjects::radiationDose::fields() const
{
    return wordList({UName_, GName_});
}


bool Foam::functionObjects::radiationDose::execute()
{
    // v0.1 skeleton: integration loop is implemented in subsequent commits
    return true;
}


bool Foam::functionObjects::radiationDose::write()
{
    // v0.1 skeleton: output writers are implemented in subsequent commits
    return true;
}


// ************************************************************************* //
