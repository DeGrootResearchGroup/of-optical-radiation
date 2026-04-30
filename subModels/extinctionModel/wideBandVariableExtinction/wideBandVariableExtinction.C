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

#include "wideBandVariableExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(wideBandVariableExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            wideBandVariableExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::wideBandVariableExtinction::wideBandVariableExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_((dict.subDict(typeName + "Coeffs"))),
    nBands_(readLabel(coeffsDict_.lookup("nBands"))),
    absorption_(readBool(coeffsDict_.lookup("absorption"))),
    nAbsorbing_(readLabel(coeffsDict_.lookup("nAbsorbing"))),
    aSpecies_(nAbsorbing_),
    aFields_(nAbsorbing_),
    aCoeffs_(nAbsorbing_),
    scattering_(readBool(coeffsDict_.lookup("scattering"))),
    nScattering_(readLabel(coeffsDict_.lookup("nScattering"))),
    sSpecies_(nScattering_),
    sFields_(nScattering_),
    sCoeffs_(nScattering_)
{
    // Initialize the extinction model
    init(nBands_);

    // Initialize absorption coefficients
    forAll(aCoeffs_, i)
    {
        aCoeffs_[i].setSize(nBands_);
        forAll(aCoeffs_[i], j)
        {
            aCoeffs_[i][j] = 0.0;
        }
    }

    // Initialize scattering coefficients
    forAll(sCoeffs_, i)
    {
        sCoeffs_[i].setSize(nBands_);
        forAll(sCoeffs_[i], j)
        {
            sCoeffs_[i][j] = 0.0;
        }
    }

    // Read the coefficients
    if (absorption_) coeffsDict_.lookup("absorptionCoeff") >> aCoeffs_;
    if (scattering_) coeffsDict_.lookup("scatteringCoeff") >> sCoeffs_;

    // Read the species names
    if (absorption_) coeffsDict_.lookup("absorbingVars") >> aSpecies_;
    if (scattering_) coeffsDict_.lookup("scatteringVars") >> sSpecies_;

    // Read the absorbing species fields'
    /*
    forAll(aFields_, i)
    {
        aFields_.set
        (
           i,
           mesh.lookupObject<volScalarField>(aSpecies_[i])
        );
    }

    // Read the scattering species fields
    forAll(sFields_, i)
    {
        sFields_.set
        (
            i,
            mesh.lookupObject<volScalarField>(sSpecies_[i])
        );
    }
    */
    // Correct the extinction coefficient fields
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::photoBio::wideBandVariableExtinction::~wideBandVariableExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::photoBio::wideBandVariableExtinction::correct()
{
    
    forAll(aFields_, i)
    {
        aFields_.set
        (
           i,
           this->mesh().lookupObject<volScalarField>(aSpecies_[i])
        );
    }

    // Read the scattering species fields
    forAll(sFields_, i)
    {
        sFields_.set
        (
            i,
            this->mesh().lookupObject<volScalarField>(sSpecies_[i])
        );
    }


    // Set the absorption coefficient field
  forAll(ALambda_, iBand)
  {
      ALambda_[iBand] = dimensionedScalar("A", dimless/dimLength, 0.0);
      for (label i = 0; i < nAbsorbing_; i++)
      {
          dimensionedScalar a("a", dimLength*dimLength/dimMass, aCoeffs_[i][iBand]);
          ALambda_[iBand] += a*aFields_[i];
      }
  }

  // Set the scattering coefficient field
  forAll(SLambda_, iBand)
  {
      SLambda_[iBand] = dimensionedScalar("S", dimless/dimLength, 0.0);
      for (label i = 0; i < nScattering_; i++)
      {
          dimensionedScalar s("s", dimLength*dimLength/dimMass, sCoeffs_[i][iBand]);
          SLambda_[iBand] += s*sFields_[i];
      }
  }
}


// ************************************************************************* //


