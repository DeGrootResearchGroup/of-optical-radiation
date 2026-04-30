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
#include "phaseFunctionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::optical::phaseFunctionModel> Foam::optical::phaseFunctionModel::New
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
{
    word phaseFunctionModelType(dict.lookup("phaseFunctionModel"));
    
    Info<< "Selecting phaseFunctionModel " << phaseFunctionModelType << endl;
    
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(phaseFunctionModelType);
    
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "phaseFunctionModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown phaseFunctionModelType type "
            << phaseFunctionModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid phaseFunctionModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

        return autoPtr<phaseFunctionModel>(cstrIter()(dom,dict,nDim));
}


// ************************************************************************* //
