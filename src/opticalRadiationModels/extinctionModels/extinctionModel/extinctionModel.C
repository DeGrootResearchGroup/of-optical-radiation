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

#include "extinctionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(extinctionModel, 0);
        defineRunTimeSelectionTable(extinctionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::extinctionModel::extinctionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::optical::extinctionModel::~extinctionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::extinctionModel::init(const label nBands)
{
  // Set the number of bands and size of pointer lists
  nBands_ = nBands;
  ALambda_.setSize(nBands_);
  SLambda_.setSize(nBands_);

  // Create absorption coefficient fields
  forAll(ALambda_, iBand)
  {
      ALambda_.set
      (
          iBand,
          new volScalarField
          (
              IOobject
              (
                  "ALambda_" + Foam::name(iBand) ,
                  mesh_.time().name(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              mesh_,
              dimless/dimLength
          )
      );
  }

  // Create scattering coefficient fields
  forAll(SLambda_, iBand)
  {
      SLambda_.set
      (
          iBand,
          new volScalarField
          (
              IOobject
              (
                  "SLambda_" + Foam::name(iBand) ,
                  mesh_.time().name(),
                  mesh_,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE
              ),
              mesh_,
              dimless/dimLength
          )
      );
  }
}


// ************************************************************************* //
