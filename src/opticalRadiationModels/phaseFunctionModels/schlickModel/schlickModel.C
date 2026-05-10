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
#include "DOM.H"

using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(schlickModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            schlickModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::schlickModel::schlickModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom,dict,nDim),
    coeffsDict_(dict.subDict(typeName + "Coeffs"))
{
    const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    functionDicts.lookup("inScatter") >> inScatter_;

    if(inScatter_)
    {
        nBand_ = dom_.nBand();
        nAngle_ = dom_.nAngle();
        functionDicts.lookup("subAngleNum") >> subAngleNum_;
        k_.setSize(nBand_);
        functionDicts.lookup("asymmetryFactor") >> k_;
        forAll(k_, b)
        {
            // Schlick (1+k^2)/(4 pi (1+k cosV)^2) goes singular at
            // 1 + k cosV = 0; for cosV in [-1, 1] this is reached at
            // |k| = 1 (k = -1 with forward cosV = 1, or k = +1 with
            // backward cosV = -1). |k| > 1 is non-physical anyway
            // (Schlick's k plays the same mean-cosine role as HG's g).
            if (mag(k_[b]) >= 1.0)
            {
                FatalErrorInFunction
                    << "asymmetryFactor[" << b << "] = " << k_[b]
                    << " is out of range; Schlick requires |k| < 1"
                    << exit(FatalError);
            }
        }
        const label nPhi = dom_.nPhi();
        const label nTheta = dom_.nTheta();   
        const scalar deltaPhi   =  pi / (2.0*nPhi);
        const scalar deltaTheta =  pi / nTheta;
        
        phaseFunction_.setSize(nAngle_*nAngle_*nBand_);    
        for(label i =0; i< nAngle_*nAngle_*nBand_;i++)     phaseFunction_[i] = 0.0;
        
        if (nDim == 3)    //3D 
        {     
            scalar dp = deltaPhi/subAngleNum_;
            scalar dt = deltaTheta/subAngleNum_;
	
            for(label iband = 0; iband < nBand_ ; iband++)
            {
		for(scalar i = 0; i<nAngle_ ; i++)
		{
                    label rayI = i + iband*nAngle_;  
                    scalar pfSum = 0;
                    for(scalar j = 0; j<nAngle_ ; j++)
                    {
                        label rayJ = j + iband*nAngle_;  
                        label idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
                        for(scalar m = 0; m < subAngleNum_ ; m++)
                        {
                            for(scalar n = 0; n < subAngleNum_ ; n++)
                            {
                                scalar nP = (2.0*m-subAngleNum_+1.0)*dp/2.0+dom.IRay(rayJ).phi();
                                scalar nT = (2.0*n-subAngleNum_+1.0)*dt/2.0+dom.IRay(rayJ).theta();
                                scalar nOmega = 2*sin(nT)*sin(dt/2)*dp;
                                vector nD = vector (sin(nT)*cos(nP), sin(nT)*sin(nP), cos(nT));
                                scalar cosV = dom.IRay(rayI).d()  & nD;

                                phaseFunction_[idx]=phaseFunction_[idx]+sl3d(cosV,k_[iband])*nOmega;
                            }
                        }
                        pfSum = pfSum + phaseFunction_(idx);
                        phaseFunction_[idx] = phaseFunction_[idx]/dom.IRay(rayJ).omega() ;
                    }
			
                    for(scalar j = 0; j<nAngle_ ; j++)
                    {
			label  idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
			phaseFunction_[idx] = phaseFunction_[idx]/pfSum;
                    }
		}
            }
	}
        
        if (nDim == 2)    //2D 
        {
            scalar  dp = deltaPhi/subAngleNum_;
            
            for(label iband = 0; iband < nBand_ ; iband++)
            {
		for(label i = 0; i<nAngle_ ; i++)
		{
                    label rayI = i + iband*nAngle_;  
                    scalar pfSum = 0;
                    for(label j = 0; j<nAngle_; j++)
                    {
                        label rayJ = j + iband*nAngle_;  
                        label idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
                        for(label m = 0; m < subAngleNum_ ; m++)
                        {
                            scalar nP = (2.0*m -subAngleNum_ +1.0)*dp/2.0 + dom.IRay(rayJ).phi(); 
                            scalar nOmega = 2*dp;
                            vector nD = vector (cos(nP), sin(nP), 0);
                            scalar cosV = dom.IRay(rayI).d() & nD;
                            phaseFunction_(idx) = phaseFunction_(idx) + sl2d(cosV,k_[iband])*nOmega;
                        }
			
                        pfSum = pfSum + phaseFunction_(idx);
                        phaseFunction_[idx] = phaseFunction_[idx]/dom.IRay(rayJ).omega() ;
                    }
                    
                    for(scalar j = 0; j<nAngle_ ; j++)
                    {
			label  idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
			phaseFunction_[idx] = phaseFunction_[idx]/pfSum;
                    }
		}
            }
	}
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::schlickModel::~schlickModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar  Foam::optical::schlickModel::correct
(
    const label rayI,
    const label rayJ,
    const label iBand
) const
{
    // See HenyeyGreensteinModel::correct() for the band-offset rationale.
    const label i = rayI - iBand*nAngle_;
    const label j = rayJ - iBand*nAngle_;
    return phaseFunction_[j + i*nAngle_ + iBand*nAngle_*nAngle_];
}

 Foam::scalar  Foam::optical::schlickModel::sl3d
(
    const scalar cosV,
    const scalar k
) const 
{
    return (1 + pow(k,2)) / (4.0 * pi * pow(1 + k*cosV, 2.0));
}


Foam::scalar  Foam::optical::schlickModel::sl2d
(
    const scalar cosV,
    const scalar k
) const
{
    return (1 + pow(k,2)) / (2.0 * pi * pow(1 + k*cosV, 2.0));

}
// ************************************************************************* //
