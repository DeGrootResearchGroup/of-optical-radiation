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
    Foam::optical::isotropicModel

Description
    Isotropic phase function model: scattering is enabled and the phase
    function P(theta) is uniformly 1 for every ray pair. Equivalent to
    HenyeyGreensteinModel with g = 0, but without the per-call pow / mul /
    sub / div arithmetic, and with a clearer dictionary intent
    ("phaseFunction { type isotropic; }").

SourceFiles
    isotropicModel.C

\*---------------------------------------------------------------------------*/

#ifndef isotropicModel_H
#define isotropicModel_H

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

        //- Construct from components
    isotropicModel
    (
        const DOM& dom,
        const dictionary& dict,
        const label& nDim
    );


    // Destructor

    virtual ~isotropicModel();


    // Member Functions

        //- Scattering is enabled (overrides base class default of false)
        virtual bool inScatter() const
        {
            return true;
        }

        //- Phase function P(theta) = 1 for all ray pairs (isotropic)
        virtual scalar correct
        (
            const label rayI,
            const label rayJ,
            const label iBand
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
