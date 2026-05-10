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
    Foam::optical::rayleighExtinction

Description
    Rayleigh scattering of an ideal gas. Per-band scattering coefficient

        sigma_s(lambda) = N(T,p) * 24 pi^3 / (lambda^4 N_s^2)
                                 * ((n^2 - 1)/(n^2 + 2))^2 * F_K

    with N(T,p) = p/(k_B T), n(lambda) the gas refractive index at standard
    conditions (Peck & Reeder 1972 fit for air, used here below 230 nm by
    extrapolation), and F_K the King correction factor (typical 1.05 for
    air). Absorption is zero -- pure scattering. Compose with another
    extinction model via "composite" if the gas absorbs (e.g. O2/O3 in the
    UV).

SourceFiles
    rayleighExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef optical_rayleighExtinction_H
#define optical_rayleighExtinction_H

#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                  Class rayleighExtinction Declaration
\*---------------------------------------------------------------------------*/

class rayleighExtinction
:
    public extinctionModel
{

    // Private data

        //- Extinction model dictionary
        dictionary coeffsDict_;

        //- Per-band central wavelengths [nm]
        List<scalar> wavelengths_;

        //- Gas temperature [K]
        scalar T_;

        //- Gas pressure [Pa]
        scalar p_;

        //- King correction factor (depolarisation), ~1.05 for air
        scalar kingFactor_;

        //- Resolved scattering coefficient per band [1/m], cached for write()
        List<scalar> sigmaSPerBand_;

public:

    //- Runtime type information
    TypeName("rayleigh");


    // Constructors

        //- Construct from components
        rayleighExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~rayleighExtinction();


    // Member Functions

        //- No-op: coefficients are uniform and set in the constructor.
        virtual void correct() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
