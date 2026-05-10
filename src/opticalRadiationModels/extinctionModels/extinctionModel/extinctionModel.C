/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

  // When constructed as a child of compositeExtinction, the per-band fields
  // are intermediate: registering them would collide with sibling and
  // composite ALambda_<i> / SLambda_<i> names, and writing them would clutter
  // time directories with quantities the composite already exposes.
  const bool isChild =
      dict_.lookupOrDefault<bool>("_compositeChild", false);
  const IOobject::writeOption wOpt =
      isChild ? IOobject::NO_WRITE : IOobject::AUTO_WRITE;
  const bool registerField = !isChild;

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
                  wOpt,
                  registerField
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
                  wOpt,
                  registerField
              ),
              mesh_,
              dimless/dimLength
          )
      );
  }
}


// ************************************************************************* //
