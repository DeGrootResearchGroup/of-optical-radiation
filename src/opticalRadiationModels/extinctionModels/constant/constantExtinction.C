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

#include "constantExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(constantExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            constantExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::constantExtinction::constantExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    absorption_(coeffsDict_.lookup("absorption")),
    scattering_(coeffsDict_.lookup("scattering")),
    ABand_(),
    SBand_()
{
    // init() sets the inherited nBands_ and resizes ALambda_ / SLambda_.
    init(readLabel(coeffsDict_.lookup("nBands")));
    ABand_.setSize(nBands_, 0.0);
    SBand_.setSize(nBands_, 0.0);

    // Read the coefficients
    if (absorption_) coeffsDict_.lookup("absorptionCoeff") >> ABand_;
    if (scattering_) coeffsDict_.lookup("scatteringCoeff") >> SBand_;

    // Set the absorption coefficient field
    forAll(ALambda_, iBand)
    {
        ALambda_[iBand] = dimensionedScalar("A", dimless/dimLength, ABand_(iBand));
    }

    // Set the scattering coefficient field
    forAll(SLambda_, iBand)
    {
        SLambda_[iBand] = dimensionedScalar("S", dimless/dimLength, SBand_(iBand));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::optical::constantExtinction::~constantExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::optical::constantExtinction::correct()
{
    // Nothing to be done for constant coefficients
    return;
}


// ************************************************************************* //


