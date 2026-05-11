/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(dispersionModel, 0);
    defineRunTimeSelectionTable(dispersionModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::dispersionModel::dispersionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::dose::dispersionModel::~dispersionModel()
{}


// * * * * * * * * * * * Static requiredFields Dispatch * * * * * * * * * * //

namespace
{
    // File-scope table for the per-type static requiredFields fn
    // pointers. Wrapped in a function to guarantee construct-on-first-
    // use ordering across translation units.
    Foam::HashTable<Foam::dose::dispersionModel::requiredFieldsFn>&
    requiredFieldsTable()
    {
        static Foam::HashTable
        <Foam::dose::dispersionModel::requiredFieldsFn> tbl;
        return tbl;
    }
}


Foam::dose::dispersionModel::addRequiredFields::addRequiredFields
(
    const word& type,
    requiredFieldsFn fn
)
{
    requiredFieldsTable().insert(type, fn);
}


Foam::wordList Foam::dose::dispersionModel::requiredFields
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
            << "No requiredFields entry registered for dispersionModel"
            << " type '" << t << "'. Registered types: "
            << tbl.sortedToc()
            << exit(FatalError);
    }
    return (*iter)(dict);
}


// ************************************************************************* //
