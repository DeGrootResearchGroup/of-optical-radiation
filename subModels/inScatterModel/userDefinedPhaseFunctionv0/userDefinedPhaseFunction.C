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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
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

Foam::optical::userDefinedPhaseFunction::userDefinedPhaseFunction
(
    const dictionary& dict
  //  , const fvMesh& mesh
)
:
 //   inScatterModel(dict, mesh),
     inScatterModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs"))
 //  , totalWaveLength_(0)
{
    label iBand = 0;
    const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    
    functionDicts.lookup("inScatter") >> inScatter_;
    functionDicts.lookup("nBand") >> nBands_;
    functionDicts.lookup("phaseFunctionAngleNum") >> phaseFunctionAngleNum_;
    
    phaseFunc_.setSize(nBand*phaseFunctionAngleNum_);   
   
    if(inScatter_)
    {
		functionDicts.lookup("phaseFunction") >> phaseFunc_;
	//	forAll(phaseFunc_,i)  phaseFunc_[i] = phaseFunc_[i]*bioDensity_;
	}else
	{
		forAll(phaseFunc_,i)  phaseFunc_[i] = 0.0;
	}
	    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::userDefinedPhaseFunction::~userDefinedPhaseFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


 Foam::scalar  Foam::optical::userDefinedPhaseFunction::correct
(
          const  scalar  rayCos,
          const   label  iBand
) const
{
/*	   scalar deltaAngle_ = 180/scalar(phaseFunctionAngleNum_-1);
       label  i0 = label(rayCos/deltaAngle_);
       scalar x0 = scalar(i0)*deltaAngle_; 
       scalar y = phaseFunc_[iBand][i0]+(phaseFunc_[iBand][i0+1]-phaseFunc_[iBand][i0])/deltaAngle_*(rayCos-x0); */
    
      label i = label(foam::acos (rayCos)/(180.0/phaseFunctionAngleNum_));   
	  return phaseFunc_[iBand*phaseFunctionAngleNum_+i];

}





// ************************************************************************* //
