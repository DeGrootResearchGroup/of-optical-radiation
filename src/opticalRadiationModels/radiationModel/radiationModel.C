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

#include "radiationModel.H"
#include "extinctionModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(radiationModel, 0);
        defineRunTimeSelectionTable(radiationModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::radiationModel::radiationModel(const volScalarField& I)
:
    IOdictionary
    (
        IOobject
        (
            "opticalRadiationProperties",
            I.time().constant(),
            I.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(I.mesh()),
    time_(I.time()),
    opticalRadiation_(false),
    coeffs_(dictionary::null),
    extinction_(nullptr)
{}


Foam::optical::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& I
)
:
    IOdictionary
    (
        IOobject
        (
            "opticalRadiationProperties",
            I.time().constant(),
            I.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(I.mesh()),
    time_(I.time()),
    opticalRadiation_(lookup("opticalRadiation")),
    coeffs_(subDict(type + "Coeffs")),
    extinction_(extinctionModel::New(*this, mesh_))
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::optical::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::optical::radiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("opticalRadiation") >> opticalRadiation_;
        coeffs_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::optical::radiationModel::correct()
{
   if (!opticalRadiation_)
    {
        return;
    }
    calculate();
}




// ************************************************************************* //
