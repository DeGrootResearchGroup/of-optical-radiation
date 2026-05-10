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

Class
    Foam::optical::HenyeyGreensteinModel

Description
    Henyey-Greenstein anisotropic phase function:

        Phi(cos theta) = (1 - g^2) / (4 pi (1 + g^2 - 2 g cos theta)^1.5)

    parameterised by a single asymmetry factor g per band (g = 0 is
    isotropic; g -> 1 is fully forward-peaked; g -> -1 is fully
    backward-peaked). Requires |g| < 1 for non-singularity.

SourceFiles
    HenyeyGreensteinModel.C

\*---------------------------------------------------------------------------*/

#ifndef HenyeyGreensteinModel_H
#define HenyeyGreensteinModel_H

#include "phaseFunctionModel.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

class DOM;

/*---------------------------------------------------------------------------*\
                    Class HenyeyGreensteinModel Declaration
\*---------------------------------------------------------------------------*/

class HenyeyGreensteinModel
:
    public phaseFunctionModel
{
private:

    // Private data

        //- Asymmetry factor g per band, |g| < 1.
        DynamicList<scalar, 1, 0> g_;


    // Private member functions

        //- Override of phaseFunctionModel::phaseShape: returns the HG
        //  shape (1 + g^2 - 2 g cos theta)^(-1.5). The (1 - g^2)/(4 pi)
        //  prefactor is omitted because the row-normalisation in
        //  buildPhaseTable absorbs it.
        virtual scalar phaseShape
        (
            const scalar cosV,
            const label iBand
        ) const override;


public:

    //- Runtime type information
    TypeName("HenyeyGreensteinModel");


    // Constructors

        HenyeyGreensteinModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~HenyeyGreensteinModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
