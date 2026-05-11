/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dose::dispersionModel>
Foam::dose::dispersionModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word modelType(dict.lookup("type"));

    Info<< "radiationDose: selecting dispersion model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown dispersionModel type " << modelType << nl << nl
            << "Valid dispersionModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<dispersionModel>(cstrIter()(dict, mesh));
}


// ************************************************************************* //
