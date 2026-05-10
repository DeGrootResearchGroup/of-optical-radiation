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
    Foam::optical::schlickModel

Description
    Schlick approximation to the Henyey-Greenstein phase function:

        Phi(cos theta) ~ 1 / (1 + k cos theta)^2

    parameterised by a single asymmetry factor k per band (k = 0 is
    isotropic; k -> 1 is forward-peaked; k -> -1 is backward-peaked).
    Requires |k| < 1 for non-singularity.

SourceFiles
    schlickModel.C

\*---------------------------------------------------------------------------*/

#ifndef optical_schlickModel_H
#define optical_schlickModel_H

#include "phaseFunctionModel.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

class DOM;

/*---------------------------------------------------------------------------*\
                          Class schlickModel Declaration
\*---------------------------------------------------------------------------*/

class schlickModel
:
    public phaseFunctionModel
{
private:

    // Private data

        //- Asymmetry factor k per band, |k| < 1.
        DynamicList<scalar, 1, 0> k_;


    // Private member functions

        //- Override of phaseFunctionModel::phaseShape: returns the Schlick
        //  shape (1 + k cos theta)^(-2). The (1 - k^2)/(4 pi) prefactor
        //  is omitted because the row-normalisation in buildPhaseTable
        //  absorbs it.
        virtual scalar phaseShape
        (
            const scalar cosV,
            const label iBand
        ) const override;


public:

    //- Runtime type information
    TypeName("schlickModel");


    // Constructors

        schlickModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~schlickModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
