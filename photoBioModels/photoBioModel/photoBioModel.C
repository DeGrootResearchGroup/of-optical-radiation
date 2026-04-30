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

#include "photoBioModel.H"
#include "extinctionModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(photoBioModel, 0);
        defineRunTimeSelectionTable(photoBioModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::photoBioModel::photoBioModel(const volScalarField& intensity)
:
    IOdictionary
    (
        IOobject
        (
            "photoBioProperties",
            intensity.time().constant(),
            intensity.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(intensity.mesh()),
    time_(intensity.time()),
    photoBio_(false),
    coeffs_(dictionary::null),
    extinction_(nullptr)
{}


Foam::photoBio::photoBioModel::photoBioModel
(
    const word& type,
    const volScalarField& intensity
)
:
    IOdictionary
    (
        IOobject
        (
            "photoBioProperties",
            intensity.time().constant(),
            intensity.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(intensity.mesh()),
    time_(intensity.time()),
    photoBio_(lookup("photoBio")),
    coeffs_(subDict(type + "Coeffs")),
    extinction_(extinctionModel::New(*this, mesh_))
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::photoBio::photoBioModel::~photoBioModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::photoBio::photoBioModel::read()
{
    if (regIOobject::read())
    {
        lookup("photoBio") >> photoBio_;
        coeffs_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::photoBio::photoBioModel::correct()
{
   if (!photoBio_)
    {
        return;
    }
    calculate();
}




// ************************************************************************* //
