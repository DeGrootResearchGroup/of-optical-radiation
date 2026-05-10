/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(opticalRadiation, 0);
    addToRunTimeSelectionTable(solver, opticalRadiation, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::opticalRadiation::opticalRadiation(fvMesh& mesh)
:
    solver(mesh),
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

Foam::solvers::opticalRadiation::~opticalRadiation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::opticalRadiation::maxDeltaT() const
{
    return GREAT;
}


void Foam::solvers::opticalRadiation::preSolve()
{
    if (mesh().changing())
    {
        FatalErrorInFunction
            << "opticalRadiation solver module does not support "
            << "dynamic meshes (mesh motion or topology change). "
            << "See the 'Mesh-motion limitations' section in "
            << "CLAUDE.md for details and the planned fix path."
            << exit(FatalError);
    }
    radiationModel_->correct();
}


bool Foam::solvers::opticalRadiation::writeData(Ostream& os) const
{
    return os.good();
}


bool Foam::solvers::opticalRadiation::read()
{
    if (solver::read())
    {
        return radiationModel_->read();
    }
    return false;
}


// ************************************************************************* //
