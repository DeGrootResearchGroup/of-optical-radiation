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
    Foam::optical::rayleighModel

Description
    Rayleigh phase function:

        Phi(cos theta) ~ 1 + cos^2 theta

    Wavelength- and band-independent. The (3/(16 pi)) prefactor of the
    canonical Rayleigh form is omitted because the row-normalisation in
    buildPhaseTable absorbs it.

SourceFiles
    rayleighModel.C

\*---------------------------------------------------------------------------*/

#ifndef optical_rayleighModel_H
#define optical_rayleighModel_H

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

class rayleighModel
:
    public phaseFunctionModel
{
private:

    // Private member functions

        //- Override of phaseFunctionModel::phaseShape: returns the
        //  Rayleigh shape (1 + cos^2 theta).
        virtual scalar phaseShape
        (
            const scalar cosV,
            const label iBand
        ) const override;


public:

    //- Runtime type information
    TypeName("rayleighModel");


    // Constructors

        rayleighModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~rayleighModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
