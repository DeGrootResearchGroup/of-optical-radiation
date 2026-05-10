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

#include "schlickModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(schlickModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            schlickModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::schlickModel::schlickModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom, dict, nDim),
    k_()
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
    k_.setSize(nBand_);
    coeffs.lookup("asymmetryFactor") >> k_;
    forAll(k_, b)
    {
        // Schlick (1 + k cosV)^2 goes singular at 1 + k cosV = 0, which
        // is reached for cosV in [-1, 1] when |k| = 1. |k| > 1 is
        // non-physical anyway (k plays the same mean-cosine role as
        // HG's g).
        if (mag(k_[b]) >= 1.0)
        {
            FatalErrorInFunction
                << "asymmetryFactor[" << b << "] = " << k_[b]
                << " is out of range; Schlick requires |k| < 1"
                << exit(FatalError);
        }
    }

    buildPhaseTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::schlickModel::~schlickModel()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::optical::schlickModel::phaseShape
(
    const scalar cosV,
    const label iBand
) const
{
    const scalar k = k_[iBand];
    const scalar denom = 1.0 + k*cosV;
    return 1.0/(denom*denom);
}


// ************************************************************************* //
