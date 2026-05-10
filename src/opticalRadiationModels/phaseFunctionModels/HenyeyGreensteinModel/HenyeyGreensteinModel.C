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

#include "HenyeyGreensteinModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(HenyeyGreensteinModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            HenyeyGreensteinModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::HenyeyGreensteinModel::HenyeyGreensteinModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom, dict, nDim),
    g_()
{
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");
    coeffs.lookup("inScatter") >> inScatter_;

    if (!inScatter_)
    {
        return;
    }

    nBand_ = dom_.nBand();
    nAngle_ = dom_.nAngle();
    coeffs.lookup("subAngleNum") >> subAngleNum_;
    g_.setSize(nBand_);
    coeffs.lookup("asymmetryFactor") >> g_;
    forAll(g_, b)
    {
        // |g| = 1 makes the denominator (1 + g^2 - 2 g cos theta)^1.5
        // singular at cos theta = 1 (forward) or -1 (backward); |g| > 1
        // is non-physical (g is the mean cosine of the scattering angle).
        if (mag(g_[b]) >= 1.0)
        {
            FatalErrorInFunction
                << "asymmetryFactor[" << b << "] = " << g_[b]
                << " is out of range; HG requires |g| < 1"
                << exit(FatalError);
        }
    }

    buildPhaseTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::HenyeyGreensteinModel::~HenyeyGreensteinModel()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::optical::HenyeyGreensteinModel::phaseShape
(
    const scalar cosV,
    const label iBand
) const
{
    const scalar g = g_[iBand];
    return 1.0/Foam::pow(1.0 + g*g - 2.0*g*cosV, 1.5);
}


// ************************************************************************* //
