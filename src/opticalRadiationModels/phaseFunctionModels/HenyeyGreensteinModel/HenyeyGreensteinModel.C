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

#include "HenyeyGreensteinModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"

using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(HenyeyGreensteinModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            HenyeyGreensteinModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::HenyeyGreensteinModel::HenyeyGreensteinModel
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
        g_.setSize(nBand_);
        functionDicts.lookup("asymmetryFactor") >> g_;
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
		for(label i = 0; i<nAngle_ ; i++)
		{
                    label rayI = i + iband*nAngle_;
                    scalar pfSum = 0;
                    for(label j = 0; j<nAngle_ ; j++)
                    {
                        label rayJ = j + iband*nAngle_;
                        label idx = j + i*nAngle_ +iband*nAngle_*nAngle_ ;
                        for(label m = 0; m < subAngleNum_ ; m++)
                        {
                            for(label n = 0; n < subAngleNum_ ; n++)
                            {
                                scalar nP = (2.0*m-subAngleNum_+1.0)*dp/2.0+dom.IRay(rayJ).phi();
                                scalar nT = (2.0*n-subAngleNum_+1.0)*dt/2.0+dom.IRay(rayJ).theta();
                                scalar nOmega = 2*sin(nT)*sin(dt/2)*dp;
                                vector nD = vector (sin(nT)*cos(nP), sin(nT)*sin(nP), cos(nT));
                                scalar cosV = dom.IRay(rayI).d()  & nD;

                                phaseFunction_[idx]=phaseFunction_[idx]+hg3d(cosV,g_[iband])*nOmega;
                            }
                        }
                        pfSum = pfSum + phaseFunction_[idx];
                        phaseFunction_[idx] = phaseFunction_[idx]/dom.IRay(rayI).omega() ;
                    }

                    for(label j = 0; j<nAngle_ ; j++)
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
                            phaseFunction_(idx) = phaseFunction_(idx) + hg2d(cosV,g_[iband])*nOmega;
                        }
			
                        pfSum = pfSum + phaseFunction_(idx);
                        phaseFunction_[idx] = phaseFunction_[idx]/dom.IRay(rayI).omega() ;
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

Foam::optical::HenyeyGreensteinModel::~HenyeyGreensteinModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar  Foam::optical::HenyeyGreensteinModel::correct
(
    const label rayI,
    const label rayJ,
    const label iBand
) const
{
    // rayI / rayJ are flat ray indices that already include the band
    // offset (rayI = iAngleI + iBand*nAngle_); the precomputed phase
    // function table is indexed by per-band angle indices, so subtract
    // the band offset before forming the index. (Without this, only
    // band 0 indexed correctly and bands >=1 read out-of-bounds garbage,
    // producing NaN scatter sources.)
    const label i = rayI - iBand*nAngle_;
    const label j = rayJ - iBand*nAngle_;
    return phaseFunction_[j + i*nAngle_ + iBand*nAngle_*nAngle_];
}

 Foam::scalar  Foam::optical::HenyeyGreensteinModel::hg3d
(
    const scalar cosV,
    const scalar g
) const 
{
// this is just fucking wrong
    //    return Foam::pow((1-g*g)/(1+g*g-2*g*cosV),1.5);
    return (1 - pow(g,2)) / (4.0 * pi * pow(1 + pow(g,2) - 2*g*cosV,1.5));
}


Foam::scalar  Foam::optical::HenyeyGreensteinModel::hg2d
(
    const scalar cosV,
    const scalar g
) const
{
    return (1 - pow(g,2)) / (2.0 * pi * pow(1 + pow(g,2) - 2*g*cosV,1.5));
}
// ************************************************************************* //
