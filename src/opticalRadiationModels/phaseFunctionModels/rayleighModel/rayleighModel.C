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

#include "rayleighModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(rayleighModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            rayleighModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::rayleighModel::rayleighModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom, dict, nDim)
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

    buildPhaseTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::rayleighModel::~rayleighModel()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::optical::rayleighModel::phaseShape
(
    const scalar cosV,
    const label /*iBand*/
) const
{
    return 1.0 + cosV*cosV;
}


// ************************************************************************* //
