/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "compositeExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(compositeExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            compositeExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::compositeExtinction::compositeExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    nBands_(readLabel(coeffsDict_.lookup("nBands"))),
    models_()
{
    init(nBands_);

    const dictionary& modelsDict = coeffsDict_.subDict("models");
    const wordList childNames = modelsDict.toc();

    if (childNames.empty())
    {
        FatalErrorIn
        (
            "compositeExtinction::compositeExtinction"
            "(const dictionary&, const fvMesh&)"
        )   << "compositeCoeffs.models contains no child models"
            << exit(FatalError);
    }

    models_.setSize(childNames.size());

    forAll(childNames, i)
    {
        Info<< "    composite child " << i << ": " << childNames[i] << endl;

        // Copy the child sub-dict so we can inject the composite-child flag
        // (suppresses field registration / AUTO_WRITE in init()).
        dictionary childDict(modelsDict.subDict(childNames[i]));
        childDict.set("_compositeChild", true);

        models_.set(i, extinctionModel::New(childDict, mesh).ptr());

        if (models_[i].nBands() != nBands_)
        {
            FatalErrorIn
            (
                "compositeExtinction::compositeExtinction"
                "(const dictionary&, const fvMesh&)"
            )   << "Child model '" << childNames[i] << "' has nBands = "
                << models_[i].nBands() << " but composite has nBands = "
                << nBands_
                << exit(FatalError);
        }
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::compositeExtinction::~compositeExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::compositeExtinction::correct()
{
    forAll(models_, i)
    {
        models_[i].correct();
    }

    forAll(ALambda_, iBand)
    {
        ALambda_[iBand] = dimensionedScalar("A", dimless/dimLength, 0.0);
        SLambda_[iBand] = dimensionedScalar("S", dimless/dimLength, 0.0);

        forAll(models_, i)
        {
            ALambda_[iBand] += models_[i].A(iBand);
            SLambda_[iBand] += models_[i].S(iBand);
        }
    }
}


// ************************************************************************* //
