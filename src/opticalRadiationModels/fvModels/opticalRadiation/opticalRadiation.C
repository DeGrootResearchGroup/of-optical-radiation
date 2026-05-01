/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

#include "opticalRadiation.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(opticalRadiation, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        opticalRadiation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::opticalRadiation::opticalRadiation
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    I_
    (
        new volScalarField
        (
            IOobject
            (
                "I",
                mesh.time().name(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    ),
    radiationModel_(Foam::optical::radiationModel::New(I_()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::opticalRadiation::~opticalRadiation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::opticalRadiation::addSupFields() const
{
    // No source contribution to host equations.
    return wordList();
}


bool Foam::fv::opticalRadiation::movePoints()
{
    return true;
}


void Foam::fv::opticalRadiation::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::opticalRadiation::mapMesh(const polyMeshMap&)
{}


void Foam::fv::opticalRadiation::distribute(const polyDistributionMap&)
{}


void Foam::fv::opticalRadiation::correct()
{
    radiationModel_->correct();
}


bool Foam::fv::opticalRadiation::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        return radiationModel_->read();
    }
    return false;
}


// ************************************************************************* //
