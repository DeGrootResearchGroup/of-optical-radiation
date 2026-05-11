/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "motionModel.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(motionModel, 0);
    defineRunTimeSelectionTable(motionModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::motionModel::motionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::dose::motionModel::~motionModel()
{}


// * * * * * * * * * * * Static requiredFields Dispatch * * * * * * * * * * //

namespace
{
    Foam::HashTable<Foam::dose::motionModel::requiredFieldsFn>&
    requiredFieldsTable()
    {
        static Foam::HashTable
        <Foam::dose::motionModel::requiredFieldsFn> tbl;
        return tbl;
    }
}


Foam::dose::motionModel::addRequiredFields::addRequiredFields
(
    const word& type,
    requiredFieldsFn fn
)
{
    requiredFieldsTable().insert(type, fn);
}


Foam::wordList Foam::dose::motionModel::requiredFields
(
    const dictionary& dict
)
{
    const word t(dict.lookup("type"));
    auto& tbl = requiredFieldsTable();
    auto iter = tbl.find(t);
    if (iter == tbl.end())
    {
        FatalErrorInFunction
            << "No requiredFields entry registered for motionModel"
            << " type '" << t << "'. Registered types: "
            << tbl.sortedToc()
            << exit(FatalError);
    }
    return (*iter)(dict);
}


// ************************************************************************* //
