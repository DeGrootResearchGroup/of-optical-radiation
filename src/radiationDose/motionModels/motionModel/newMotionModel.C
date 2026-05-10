/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "motionModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dose::motionModel>
Foam::dose::motionModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    const word modelType(dict.lookup("type"));

    Info<< "Selecting motion model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown motionModel type " << modelType << nl << nl
            << "Valid motionModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<motionModel>(cstrIter()(dict, mesh));
}


// ************************************************************************* //
