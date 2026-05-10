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

Class
    Foam::optical::isotropicModel

Description
    Isotropic phase function: Phi(cos theta) = constant. Bit-for-bit
    equivalent to HenyeyGreensteinModel with g = 0, but with no
    asymmetry-factor knob to set. Routes through the base class's
    buildPhaseTable so the row-normalised table value table[i,j] equals
    omega_j/(4 pi), the same shape every other phase function model
    produces in the isotropic limit.

SourceFiles
    isotropicModel.C

\*---------------------------------------------------------------------------*/

#ifndef optical_isotropicModel_H
#define optical_isotropicModel_H

#include "phaseFunctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

class DOM;

/*---------------------------------------------------------------------------*\
                       Class isotropicModel Declaration
\*---------------------------------------------------------------------------*/

class isotropicModel
:
    public phaseFunctionModel
{

public:

    //- Runtime type information
    TypeName("isotropicModel");


    // Constructors

        isotropicModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~isotropicModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
