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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::optical::inScatterModel> Foam::optical::inScatterModel::New
(
    const dictionary& dict
 //   ,const fvMesh& mesh
)
{
    word inScatterModelType(dict.lookup("inScatterModel"));

    Info<< "Selecting inScatterModel " << inScatterModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(inScatterModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "inScatterModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown inScatterModelType type "
            << inScatterModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid inScatterModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

 //   return autoPtr<inScatterModel>(cstrIter()(dict, mesh));
        return autoPtr<inScatterModel>(cstrIter()(dict));
}


// ************************************************************************* //
