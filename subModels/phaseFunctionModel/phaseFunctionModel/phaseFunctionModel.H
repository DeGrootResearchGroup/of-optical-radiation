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
    Foam::optical::scatterModel

Description
    Base class for opticalRadiation scattering

\*---------------------------------------------------------------------------*/

#ifndef phaseFunctionModel_H
#define phaseFunctionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{
class DOM;

/*---------------------------------------------------------------------------*\
                           Class phaseFunctionModel Declaration
\*---------------------------------------------------------------------------*/

class phaseFunctionModel
{
protected:
    // Protected data
    const DOM& dom_;
    
    //- opticalRadiation model dictionary
    const dictionary dict_;

    //- number of dimensions
    label nDim_;
    
    label nBand_;
    
    label nAngle_;
        
public:
    //- Runtime type information
    TypeName("phaseFunctionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
	autoPtr,
	phaseFunctionModel,
	dictionary,
	(
	    const DOM& dom,
	    const dictionary& dict,
	    const label& nDim
        ),
	(dom,dict,nDim)
    );
    
    
    // Constructors
    
    //- Construct from components
    phaseFunctionModel
    (
	const DOM& dom,
	const dictionary& dict,
	const label& nDim
    );

    // Selector
    static autoPtr<phaseFunctionModel> New
    (
	const DOM& dom,
	const dictionary& dict,
	const label& nDim
        );
    
    // Destructor
    virtual ~phaseFunctionModel();
    
    // Member Functions
    
    // Access

    //- Reference to the dictionary
    inline const dictionary& dict() const
    {
	return dict_;
    }
    
    virtual bool inScatter() const;
    
    
    virtual scalar correct
    (
	label rayI,
	label rayJ,
	label iBand
	) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

