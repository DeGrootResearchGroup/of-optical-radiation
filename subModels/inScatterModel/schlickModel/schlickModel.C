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

#include "schlickModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(schlickModel, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
	    schlickModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::schlickModel::schlickModel
(
    const dictionary& dict
  //  , const fvMesh& mesh
)
:
     inScatterModel(dict),
    coeffsDict_(dict.subDict(typeName + "Coeffs"))
{
     coeffsDict_.lookup("subAngleNum") >> subAngleNum;
     
     const photoBioModel& photoBio = db().lookupObject<photoBioModel>("photoBioProperties");
     const photoBioDOM& dom(refCast<const photoBioDOM>(photoBio));
     
     nBands_ = dom.nBands();   //    coeffsDict_.lookup("nBand") >> nBands_;
     asymmetryFactor_.setSize(nBands_);   
     coeffsDict_.lookup("asymmetryFactor") >> k_;
     
     nAngle_ = dom.nAngle();
     const  label nPhi = dom.nPhi();  
     const  label nTheta = dom.nTheta();   
     
     phaseFunction_.setSize(nAngle_*nAngle_*nBands_);    
     
     
    if (internalField().mesh().nSolutionD() == 3)    //3D 
    {
	scalar  dp = deltaPhi/subAngleNum;
	scalar  dt = deltaTheta/subAngleNum;
	
	for(iband = 0; iband < nBands_ ; iband++)
	{
		for(scalar i = 0; i<nAngle ; i++)
		{
			label rayI = i + iBand*nAngle;  
			scalar pfSum = 0;
			for(scalar j = 0; j<nAngle ; j++)
			{
				label  rayJ = j + iBand*nAngle;  
				lable  idx = j + i*nAngle +iBand*nAngle*nAngle ;
				for(scalar m = 0; m < subAngleNum ; m++)
				{
				for(scalar n = 0; n < subAngleNum ; n++)
				{
					scalar nP = (2.0*m -subAngleNum -1.0)*dp/2.0 + dom.IRay(rayJ).phi(); 
					scalar nT = (2.0*n -subAngleNum -1.0)*dt/2.0 + dom.IRay(rayJ).theta();
					scalar nOmega = 2*sin(nT)*sin(dt/2)*dp;
					vector nD = vector (sin(nT)*cos(nP), sin(nT)*sin(nP), cos(nT));
					scalar cosV = dom.IRay(rayI).d()  & nD;

					phaseFunction_(idx) = phaseFunction_(idx) + hg3d(cosV,g)*nOmega;
				}
				}
			pfSum = pfSum + phaseFunction_(idx);
			phaseFunction_(idx) = phaseFunction_(idx)/dom.IRay(rayI).omega() ;
			}
			
			for(scalar j = 0; j<nAngle ; j++)
			{
			lable  idx = j + i*nAngle +iBand*nAngle*nAngle ;
			phaseFunction_(idx) = phaseFunction_(idx)*4*pi./pfSum;
			}
		}
	}
	}
     
    if (internalField().mesh().nSolutionD() == 2)    //2D 
    {
	scalar  dp = deltaPhi/subAngleNum;
	
	for(label iband = 0; iband < nBands_ ; iband++)
	{
		for(label i = 0; i<nAngle ; i++)
		{
			label rayI = i + iBand*nAngle;  
			scalar pfSum = 0;
			for(label j = 0; j<nAngle ; j++)
			{
				label  rayJ = j + iBand*nAngle;  
				lable  idx = j + i*nAngle +iBand*nAngle*nAngle ;
				for(label m = 0; m < subAngleNum ; m++)
				{
					scalar nP = (2.0*m -subAngleNum -1.0)*dp/2.0 + dom.IRay(rayJ).phi(); 
					scalar nOmega = 2*dp;
					vector nD = vector (cos(nP), sin(nP), 0);
					scalar cosPsi = dom.IRay(rayI).d() & nD;
					phaseFunction_(idx) = phaseFunction_(idx) + hg2d(cosPsi,g(iBand))*nOmega;
				}
				}
			pfSum = pfSum + phaseFunction_(idx);
			phaseFunction_(idx) = phaseFunction_(idx)/dom.IRay(rayI).omega() ;
			}
			
			for(scalar j = 0; j<nAngle ; j++)
			{
			lable  idx = j + i*nAngle +iBand*nAngle*nAngle ;
			phaseFunction_(idx) = phaseFunction_(idx)*4*pi./pfSum;
			}
		}
	}
	}
      

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::photoBio::schlickModel::~schlickModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar  Foam::photoBio::schlickModel::correct
(
          const   label rayI,
          const   label rayJ,
          const   label  iBand
) const
{
    return phaseFunction_(rayJ + rayI*nAngle +iBand*nAngle*nAngle);
}

 Foam::scalar  Foam::photoBio::HenyeyGreensteinModel::hg3d
(
	const scalar cosPsi,
	const scalar k
) const 
{
  return (1 - pow(k,2)) / ((1 + k * cosPsi).^2)
  
    //  return (1-)/(1+g*g-2*g*cosPsi).^1.5; this is what this numpty had

}

 Foam::scalar  Foam::photoBio::HenyeyGreensteinModel::hg2d
(
	const scalar cosPsi,
	const scalar k
) const
{
  // return 0.5/pi*(1-g*g)/(1+g*g-2*g*cosV)
 return (1 - pow(k,2)) / ((1 + k * cosPsi).^2)

}



// ************************************************************************* //
