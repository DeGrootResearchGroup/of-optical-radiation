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
    FatalErrorInFunction
        << "opticalRadiation does not support moving meshes. "
        << "See the 'Mesh-motion limitations' section in CLAUDE.md "
        << "for details and the planned fix path."
        << exit(FatalError);
    return true;
}


void Foam::fv::opticalRadiation::topoChange(const polyTopoChangeMap&)
{
    FatalErrorInFunction
        << "opticalRadiation does not support topology changes "
        << "(adaptive refinement, cell add/remove). "
        << "See the 'Mesh-motion limitations' section in CLAUDE.md."
        << exit(FatalError);
}


void Foam::fv::opticalRadiation::mapMesh(const polyMeshMap&)
{
    FatalErrorInFunction
        << "opticalRadiation does not support mesh remapping "
        << "(post mesh-motion / topology-change field mapping). "
        << "See the 'Mesh-motion limitations' section in CLAUDE.md."
        << exit(FatalError);
}


void Foam::fv::opticalRadiation::distribute(const polyDistributionMap&)
{
    // Parallel redistribution (e.g. dynamic load balancing) is not
    // mesh motion -- opticalRadiation works correctly across a
    // static decomposition. No-op is intentional.
}


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
