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

Foam::autoPtr<Foam::optical::phaseFunctionModel>
Foam::optical::phaseFunctionModel::New
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
{
    // Phase function selection is optional. If the dictionary entry is
    // omitted, fall back to the base class -- inScatter_ is false by
    // default, so DOM skips the in-scatter source entirely. This is the
    // intended setup for non-scattering media.
    if (!dict.found("phaseFunctionModel"))
    {
        Info<< "No phaseFunctionModel selected; scattering disabled" << endl;
        return autoPtr<phaseFunctionModel>
        (
            new phaseFunctionModel(dom, dict, nDim)
        );
    }

    const word phaseFunctionModelType(dict.lookup("phaseFunctionModel"));

    Info<< "Selecting phaseFunctionModel " << phaseFunctionModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(phaseFunctionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseFunctionModel type "
            << phaseFunctionModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid phaseFunctionModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<phaseFunctionModel>(cstrIter()(dom, dict, nDim));
}


// ************************************************************************* //
