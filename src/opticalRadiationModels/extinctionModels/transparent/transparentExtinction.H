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
    Foam::optical::transparentExtinction

Description
    Convenience extinction model: zero absorption and scattering for all
    bands. Equivalent to "constant" with all coefficients set to zero.
    With no "transparentCoeffs" subdictionary the model is single-band;
    add "transparentCoeffs { nBands N; }" to match a multi-band DOM.

SourceFiles
    transparentExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef opticalRadiationTransparentExtinction_H
#define opticalRadiationTransparentExtinction_H

#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                  Class transparentExtinction Declaration
\*---------------------------------------------------------------------------*/

class transparentExtinction
:
    public extinctionModel
{

    public:

        //- Runtime type information
        TypeName("transparent");


    // Constructors

        //- Construct from components
        transparentExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    // Destructor
    virtual ~transparentExtinction();

    // Member functions
    virtual void correct() {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
