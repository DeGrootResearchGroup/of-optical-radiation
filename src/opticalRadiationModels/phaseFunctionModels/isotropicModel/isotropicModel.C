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

#include "isotropicModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(isotropicModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            isotropicModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::isotropicModel::isotropicModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom, dict, nDim)
{
    // No knobs: scattering is on, sub-pixel sampling is fixed at 1.
    // The base class's default phaseShape() returns 1.0 (constant), which
    // after row-normalisation gives table[i, j] = omega_j/(4 pi) -- the
    // expected isotropic shape.
    inScatter_ = true;
    nBand_ = dom_.nBand();
    nAngle_ = dom_.nAngle();
    subAngleNum_ = 1;

    buildPhaseTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::isotropicModel::~isotropicModel()
{}


// ************************************************************************* //
