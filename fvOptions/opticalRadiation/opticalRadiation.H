/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::photoBio

Description

SourceFiles
    photoBioI.H
    photoBio.C
    photoBioIO.C

\*---------------------------------------------------------------------------*/

#ifndef photoBio_H
#define photoBio_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"
#include "photoBioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
/*---------------------------------------------------------------------------*\
                         Class photoBio Declaration
\*---------------------------------------------------------------------------*/

class photoBio
:
    public option
{
    // Private data_
        //- The photoBio model pointer
        autoPtr<Foam::photoBio::photoBioModel> photoBio_;

            // Private Member Functions

        photoBio(const photoBio&) = delete;
        void operator=(const photoBio&) = delete;

public:

    //- Runtime type information
    TypeName("photoBio");

    // Constructors
        //- Construct from explicit source name and mesh
        photoBio
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Destructor
            ~photoBio();

    // IOstream Operators
    virtual bool read(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
