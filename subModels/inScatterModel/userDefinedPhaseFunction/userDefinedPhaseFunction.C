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

#include "userDefinedPhaseFunction.H"
#include "addToRunTimeSelectionTable.H"
#include "photoBioDOM.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(userDefinedPhaseFunction, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
            userDefinedPhaseFunction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::userDefinedPhaseFunction::userDefinedPhaseFunction
(
    const dictionary& dict
  //  , const fvMesh& mesh
)
:
 //   inScatterModel(dict, mesh),
    inScatterModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    nBand_(1),
    nAngle_(1),
    nDim(3)
{
    const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    
    functionDicts.readEntry("inScatter", inScatter_);

    functionDicts.readEntry("phaseFunctionAngleNum", phaseFunctionAngleNum_);
    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::photoBio::userDefinedPhaseFunction::~userDefinedPhaseFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

 Foam::scalar  Foam::photoBio::userDefinedPhaseFunction::update
(
			const	label  nPhi,
			const 	label  nTheta,
			const 	label  nBand,
			const 	label  nDim
) const
{
	nBand_ = nBand;
	nDim_  = nDim;
	nAngle_ =  4*nPhi_*nTheta_;
    phaseFunc_.setSize(nBand_*phaseFunctionAngleNum_);   
     
	if(inScatter_)
    {
		functionDicts.readEntry("phaseFunction", phaseFunc_);
	//	Info << "phaseFunc_  " << phaseFunc_ << endl;
	}
       
}



 Foam::scalar  Foam::photoBio::userDefinedPhaseFunction::correct
(
          const  scalar  rayCos,
          const   label  iBand
) const
{
	   scalar deltaAngle = Foam::constant::mathematical::pi/scalar(phaseFunctionAngleNum_-1);   
       scalar angle = 	Foam::acos(rayCos);
       label  i0 = label(angle/deltaAngle+0.2);
       
       if(i0 >= 0 && i0 < phaseFunctionAngleNum_)
       { 
               return  phaseFunc_[i0 + iBand*phaseFunctionAngleNum_];
			//	return 0.0;
		}
		else
		{ 
			return 0.0;
		}
       
}





// ************************************************************************* //
