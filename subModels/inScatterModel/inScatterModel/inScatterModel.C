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

#include "error.H"
#include "inScatterModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(inScatterModel, 0);
        defineRunTimeSelectionTable(inScatterModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::inScatterModel::inScatterModel
(
    const dictionary& dict
  //  ,const fvMesh& mesh
)
:
   //     mesh_(mesh),
        dict_(dict)
{}

// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::photoBio::inScatterModel::~inScatterModel()
{}


// ************************************************************************* //


bool  Foam::photoBio::inScatterModel::inScatter() const
{
return false;
}

void  Foam::photoBio::inScatterModel::update
(
          const  label   nPhi,
          const   label  nTheta,
          const   label  nBand,
          const   label  nDim
) const
{
}

Foam::scalar  Foam::photoBio::inScatterModel::correct
(
          const  label  rayI,
          const   label  rayJ,
          const   label  iBand
) const
{
return 0.0;
}

/*
Foam::scalar  Foam::photoBio::inScatterModel::correct
(
          const  scalar  angle,
          const   label  iBand
) const
{
return angle/10;
}
*/
