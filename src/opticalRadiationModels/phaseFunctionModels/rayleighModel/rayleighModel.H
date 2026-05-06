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
    Foam::optical::rayleighModel

Description
    Rayleigh phase function P(cos theta) ~ 1 + cos^2 theta. Wavelength- and
    band-independent. The precomputed table is row-renormalised in the same
    pixel-averaged scheme used by HenyeyGreensteinModel, so the absolute
    kernel scaling does not affect the in-scatter source -- only the angular
    shape matters.

SourceFiles
    rayleighModel.C

\*---------------------------------------------------------------------------*/

#ifndef rayleighModel_H
#define rayleighModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"

#include "phaseFunctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

class DOM;

/*---------------------------------------------------------------------------*\
                  Class rayleighModel Declaration
\*---------------------------------------------------------------------------*/

class rayleighModel : public phaseFunctionModel
{
private:

    // Private data

        Switch inScatter_;
        dictionary coeffsDict_;
        label subAngleNum_;
        DynamicList<scalar, 1, 0> phaseFunction_;

public:

    //- Runtime type information
    TypeName("rayleighModel");


    // Constructors

        //- Construct from components
        rayleighModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~rayleighModel();


    // Member Functions

        virtual bool inScatter() const
        {
            return inScatter_;
        }

        Foam::scalar correct
        (
            const label rayI,
            const label rayJ,
            const label iBand
        ) const;

        Foam::scalar rayleigh3d(const scalar cosV) const;
        Foam::scalar rayleigh2d(const scalar cosV) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
