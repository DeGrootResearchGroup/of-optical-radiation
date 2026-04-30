/*---------------------------------------------------------------------------*\
=========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "transExteriorSurfaceMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "photoBioDOM.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(0.0),
    nOwn_(0.0), 
    I0_(0.0), 
    nBands_(1),
    diffuseFraction_(0.0),
    beamWidthPhi_(0.0),
    beamWidthTheta_(0.0),
    beamDir_(vector::zero),
    beamNormToSurf_(false),
 //   photoBioBandDist_(null),
    initFlag_(0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const transExteriorSurfaceMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    diffuseFraction_(ptf.diffuseFraction_),
    beamWidthPhi_(ptf.beamWidthPhi_),
    beamWidthTheta_(ptf.beamWidthTheta_),
    beamDir_(ptf.beamDir_),
    beamNormToSurf_(ptf.beamNormToSurf_),
 //   photoBioBandDist_(ptf.photoBioBandDist_),
    initFlag_(ptf.initFlag_)
{}


Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(readScalar(dict.lookup("nNbg"))), //scalarField("n1", dict),
    nOwn_(readScalar(dict.lookup("nOwn"))), // scalarField("n2", dict);
    I0_(readScalar(dict.lookup("irradiation"))),
    nBands_(readLabel(dict.lookup("nBand"))), // scalarField("n2", dict);
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction"))),
    beamWidthPhi_(readScalar(dict.lookup("beamWidthPhi"))),
    beamWidthTheta_(readScalar(dict.lookup("beamWidthTheta"))),
    beamDir_(dict.lookup("beamDir")),
    initFlag_(0)
{
	
   photoBioBandDist_.setSize(nBands_);   

   dict.lookup("photoBioBandDist") >> photoBioBandDist_;
   dict.lookup("beamNormToSurf") >> beamNormToSurf_;

 /*     Info << "\n nNbg_  \t" << nNbg_ << endl;
   Info << "\n nOwn_  \t" << nOwn_ << endl;
    Info << "\n I0_  \t" << I0_ << endl;
    Info << "\n nBands_  \t" << nBands_ << endl;
    Info << "\n diffuseFraction_  \t" << diffuseFraction_ << endl;
    Info << "\n beamWidthPhi_  \t" << beamWidthPhi_ << endl;
    Info << "\n beamDir_  \t" << beamDir_ << endl;
    Info << "\n beamNormToSurf_  \t" << beamNormToSurf_ << endl;

    */
	    	  
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );               
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
              
    }
    else
    {
        // No value given. Restart as fixedValue b.c.
        refValue() = 0.0;
        refGrad() = 0.0;
        valueFraction() = 1.0;
        fvPatchScalarField::operator=(refValue());       
    }
    
}



Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::
transExteriorSurfaceMixedFvPatchScalarField
(
    const transExteriorSurfaceMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    diffuseFraction_(ptf.diffuseFraction_),
    beamWidthPhi_(ptf.beamWidthPhi_),
    beamWidthTheta_(ptf.beamWidthTheta_),
    beamDir_(ptf.beamDir_),
    beamNormToSurf_(ptf.beamNormToSurf_),
 //   photoBioBandDist_(ptf.photoBioBandDist_),
    initFlag_(ptf.initFlag_)
{}






// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::
updateCoeffs()
{
	if (this->updated())
    {
        return;
    }
    
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

     	
    scalarField& Iw = *this;
   
    const photoBioModel& photoBio = db().lookupObject<photoBioModel>("photoBioProperties");

    const photoBioDOM& dom(refCast<const photoBioDOM>(photoBio));
 
    if (dom.nBand() == 0)
    {
        FatalErrorIn
        (
            "Foam::photoBio::"
            "wideBandDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs"
        )   << " a non-grey boundary condition is used with a grey "
            << "absorption model" << nl << exit(FatalError);
    }
    const label patchI = patch().index();
    vectorField n = patch().Sf()/patch().magSf();
    
    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);
    
    const  label nAngle = dom.nAngle();    
    const  label iBand = dom.IRay(rayId).iBand();
    const  vector& bdRayDir = dom.IRay(rayId).d();
    const  vector& bdRayDAve = dom.IRay(rayId).dAve();
    const  scalar bdOmega = dom.IRay(rayId).omega();
    const  scalar bdRayPhi  = dom.IRay(rayId).phi();
    const  scalar bdRayTheta = dom.IRay(rayId).theta();

    const  label nTheta = dom.nTheta();    
    const  label nPhi = dom.nPhi();    
    const  scalar deltaPhi   =  pi /(2.0*nPhi);
    const  scalar deltaTheta =  pi  /nTheta;
    
    label  npPhi  = dom.nPixelPhi();
    label  npTheta = dom.nPixelTheta();
	
	label  numDivPhi =1;
	label  numDivTheta =1;
    label  rayJ = -1;
	scalar deltaDivTheta ;
	scalar beamThetaStart ;
	
	scalar nRatio = nOwn_/  nNbg_;
/*	
	
    angleDist.setSize(nPhi); 
     scalar sumP =0;      
	 for(label i = 0;i < nPhi;i++)
	 {	
		scalar alpha = i*deltaPhi ;
		angleDist[i] = Foam::cos(alpha) ;  
	    sumP = sumP + angleDist[i];

	}
	sumP = 2*sumP - angleDist[0];
	for(label i = 0;i < nPhi;i++)
	 {	angleDist[i]  = angleDist[i] /sumP;}
*/	 
	         
	 dirToAngle(beamDir_, beamAnglePhi_,beamAngleTheta_) ;
	    
     // calculate the refraction photoBio only once at the beginning
     if(initFlag_ == 0)
     {	

		  initFlag_ = 1;
		  specularRefraction.set(new scalarField(patch().size(),0.0));
		  diffusiveRefraction.set(new scalarField(patch().size(),0.0));

		 if(beamNormToSurf_ )
         {
			forAll(Iw, faceI)
			{
				scalar T = 4*nNbg_*nOwn_/(nNbg_+nOwn_)/(nNbg_+nOwn_);
				
				if(diffuseFraction_ > 0)
				{         
					diffusiveRefraction()[faceI] = diffusiveRefraction()[faceI] + I0_*T*mag(n[faceI] & bdRayDAve);  
				}
					
				dom.dirToRayId(-n[faceI], iBand, rayJ);
				if(rayJ == rayId)
				{
					specularRefraction()[faceI] = specularRefraction()[faceI] + I0_*T*bdOmega;
				}   
			} 
	      }
         else
         {

			scalar  bwP = beamWidthPhi_/180*pi;	  
			scalar  bwT = beamWidthTheta_/180*pi;
			
			if(internalField().mesh().nSolutionD() == 3)    //3D
			{		
				              
                if(bwP > deltaPhi)  
                {
					numDivPhi = npPhi*label(bwP/deltaPhi);
				}
				else
				{
					if(bwP > deltaPhi/2)  
					{
					  numDivPhi = npPhi;
					}
					else
					{
					  numDivPhi = 1;
					}
				}
				
	      
				if(bwT > deltaTheta)  
                {
					numDivTheta = npTheta*label(bwT/deltaTheta);;
				}
				else
				{
					if(bwT > deltaTheta/2)  
					{
					  numDivTheta = npTheta;
					}
					else
					{
					  numDivTheta = 1;
					}
				}
				
				 deltaDivTheta = bwT/numDivTheta ;
				 beamThetaStart = beamAngleTheta_ - bwT/2.0;;
		 }
		 else
		 {
			if (internalField().mesh().nSolutionD() == 2)    //2D (X & Y)
			{
				
			    if(bwP > deltaPhi)  
                {
					  numDivPhi = npPhi*label(bwP/deltaPhi);
				}
				else
				{
					if(bwP > deltaPhi/2)  
					{
					  numDivPhi = npPhi;
					}
					else
					{
					  numDivPhi = 1;
					}
				}
				  numDivTheta = 1;
				
				  deltaDivTheta = pi;
				  beamThetaStart = 0;
				  		
            }
			else    //1D (X)
			{
				numDivPhi = 1;	
				numDivTheta = 1;
				deltaDivTheta = pi;
				beamThetaStart = 0;
			}
		} 
		
		scalar deltaDivPhi = bwP/numDivPhi;
		scalar beamPhiStart = beamAnglePhi_ - bwP/2.0;

        forAll(Iw, faceI)
		{
			
          for (label i = 1; i <= numDivTheta; i++)
          {
			scalar  sweepDivTheta = beamThetaStart + (2*i-1)*deltaDivTheta/2.0;
          
			for (label j = 1; j <= numDivPhi; j++)
			{
				scalar  sweepDivPhi = beamPhiStart + (2*j-1)*deltaDivPhi/2.0;
          
				scalar sinTheta = Foam::sin(sweepDivTheta);
				scalar cosTheta = Foam::cos(sweepDivTheta);
				scalar sinPhi = Foam::sin(sweepDivPhi);
				scalar cosPhi = Foam::cos(sweepDivPhi);
                
					vector sweepDivDir = vector(sinTheta*cosPhi , sinTheta*sinPhi , cosTheta);
					
					vector surfNorm = -n[faceI] ;
					scalar cosA = surfNorm& sweepDivDir; 
					
					if(  cosA > 0.0 && cosA*cosA > 1-1/(nRatio*nRatio))     // direction out of the wall   
					{ 
				
					vector refracDir = (sweepDivDir - cosA*surfNorm )/nRatio
					+ Foam::sqrt(1 - 1/(nRatio*nRatio)*(1- cosA*cosA)) *surfNorm;
					
					scalar cosB = surfNorm& refracDir;  
          
					scalar r1 = (nNbg_*cosB - nOwn_*cosA )
								/  (nNbg_*cosB + nOwn_*cosA );
					scalar r2 = (nNbg_*cosA - nOwn_*cosB )
								/  (nNbg_*cosA + nOwn_*cosB );
            
					scalar T = 1.0 - 0.5*(r1*r1 + r2*r2);
          			
          			scalar refracPhi;   scalar refracTheta; 
          		    dirToAngle(refracDir, refracPhi, refracTheta) ; 
						 
					sinTheta = Foam::sin(refracTheta);
					cosTheta = Foam::cos(refracTheta);
					sinPhi = Foam::sin(refracPhi);
					cosPhi = Foam::cos(refracPhi);
					
					if(diffuseFraction_ > 0)
					{         
						vector refracdAve = vector
						(
						  cosPhi
						*Foam::sin(0.5*deltaDivPhi)
						*(deltaDivTheta - Foam::cos(2.0*refracTheta)
						*Foam::sin(deltaTheta)),     
						sinPhi
						*Foam::sin(0.5*deltaDivPhi)
						*(deltaDivTheta - Foam::cos(2.0*refracTheta)
						*Foam::sin(deltaDivTheta)),   
						0.5*deltaDivPhi
						*Foam::sin(2.0*refracTheta)*Foam::sin(deltaDivTheta)
						);
						diffusiveRefraction()[faceI] = diffusiveRefraction()[faceI] + I0_*T*mag(n[faceI] & refracdAve);  
					}
					
	                label refracRay = -1;
					dom.dirToRayId(refracDir, iBand, refracRay);
					if(refracRay == rayId)
					{

						scalar refracOmega = 2.0*sinTheta*Foam::sin(deltaDivTheta/2.0)*deltaDivPhi;			 
						specularRefraction()[faceI] = specularRefraction()[faceI] + I0_*T*refracOmega;

					}
				}
			 }
			}
		 }
		}
	}
		
	scalar specularReflection = 0.0;
	scalar diffusiveReflection =0.0;

	if (internalField().mesh().nSolutionD() == 2)    //2D (X & Y)
	{	
		npTheta = 1;
    }
    if (internalField().mesh().nSolutionD() == 1)    //2D (X & Y)
	{	
		npTheta = 1; npPhi =1;
    }
	
	     
    forAll(Iw, faceI)
    {
         specularReflection = 0.0;
         diffusiveReflection = 0.0;
	     vector surfNorm = -n[faceI] ;
	         scalar cosA = surfNorm& bdRayDir; 
	              
	   if(cosA> 0.0 )    // direction out of the wall   
       {
			if( diffuseFraction_ > 0)
			{	
				for (label jAngle = 0; jAngle < nAngle; jAngle++)
				{         
				label sweepRayID = jAngle + (iBand)*nAngle;
				vector sweepDir = dom.IRay(sweepRayID).d();
				vector sweepdAve = dom.IRay(rayJ).dAve();   
				
				scalar cosB = surfNorm& sweepDir; 
			
				if(  cosB > 0.0 )      // direction out of the wall   
				{ 
				    vector reflecIncidentDir = sweepDir - 2*cosB*surfNorm;		
					label reflecIncidentRay = -1;
                    dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                    const scalarField&  reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];     		
                     
                     if(cosB*cosB > 1 - 1/(nRatio*nRatio) )        
                     {   
						vector refracIncidentDir = (sweepDir - cosB*surfNorm )*nRatio 
							+ Foam::sqrt(1 -(nRatio*nRatio)*(1- cosB*cosB))*surfNorm ;

						scalar cosA = surfNorm& refracIncidentDir;  
					    
						scalar r1 = (nNbg_*cosB - nOwn_*cosA )
								/  (nNbg_*cosB + nOwn_*cosA );
						scalar r2 = (nNbg_*cosA - nOwn_*cosB )
								/  (nNbg_*cosA + nOwn_*cosB );
            
						scalar R = 0.5*(r1*r1 + r2*r2);
          		
						label refracIncidentRay = -1;  
						dom.dirToRayId(refracIncidentDir, 0, refracIncidentRay);	
          		    
						diffusiveReflection = diffusiveReflection + reflecFace[faceI]*R*mag(surfNorm & sweepdAve) ;  
					 }
					 else
					 {
						 diffusiveReflection = diffusiveReflection + reflecFace[faceI]*mag(surfNorm & sweepdAve) ;  
					 }

				}
				}
			}


		for (label i = 1; i <= npTheta; i++)
        {
			scalar pxRayTheta = bdRayTheta - 0.5*deltaTheta + 0.5*(2*i -1)*deltaTheta/npTheta;
				    
			for (label j = 1; j <= npPhi; j++)
			{
				scalar pxRayPhi = bdRayPhi - 0.5*deltaPhi + 0.5*(2*j -1)*deltaPhi/npPhi;
					
				scalar sinTheta = Foam::sin(pxRayTheta);
				scalar cosTheta = Foam::cos(pxRayTheta);
				scalar sinPhi = Foam::sin(pxRayPhi);
				scalar cosPhi = Foam::cos(pxRayPhi);
					
				vector  pixelDir = vector(sinTheta*cosPhi , sinTheta* sinPhi , cosTheta);
					
				scalar cosB =  pixelDir & surfNorm;
					
				if ( cosB > 0.0)   
				{
						
					scalar  pixelOmega = 2.0*sinTheta*Foam::sin(deltaTheta/2.0/npTheta)*deltaPhi/npPhi;

					vector reflecIncidentDir = pixelDir - 2*cosB*surfNorm;		
					
					label reflecIncidentRay = -1;
                    dom.dirToRayId(reflecIncidentDir, iBand, reflecIncidentRay);
                    
                    const scalarField&  reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];                                	

                    
                    if( cosB*cosB > 1-1/(nRatio*nRatio) )      // direction out of the wall   
					{ 
                    
					vector refracIncidentDir = (pixelDir - cosB * surfNorm )*nRatio 
					+ Foam::sqrt(1 -(nRatio*nRatio)*(1- cosB*cosB))*surfNorm ;
					
					scalar cosA = surfNorm& pixelDir; // / mag(d);                //    /mag(n[faceI])
          
					scalar r1 = (nNbg_*cosB - nOwn_*cosA )
								/  (nNbg_*cosB + nOwn_*cosA );
					scalar r2 = (nNbg_*cosA - nOwn_*cosB )
								/  (nNbg_*cosA + nOwn_*cosB );
            
					scalar R =  0.5*(r1*r1 + r2*r2);
          			
          			label refracIncidentRay = -1;  
          		    dom.dirToRayId(refracIncidentDir, 0, refracIncidentRay);	
 
					specularReflection = specularReflection + reflecFace[faceI]*R*pixelOmega;
					
					}
					
					else
					{
						specularReflection = specularReflection + reflecFace[faceI]*pixelOmega;   
					}
				}
			}
			}

             
       //      label  i0 = label(Foam::acos(cosA)/deltaPhi);                   
           refValue()[faceI] = diffuseFraction_*(diffusiveRefraction()[faceI] + diffusiveReflection )/pi/2     //*angleDist[i0] 
                        + (1.0 - diffuseFraction_)*(specularReflection + specularRefraction()[faceI])/bdOmega; 
                 
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            
            
        }
        else
        {
            // direction into the wall
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used
        }
    }
     
    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs(); 
    
}

void Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::dirToAngle
(
    const vector& dir,
    scalar&	tPhi,
    scalar& tTheta
) const
{

	tTheta = Foam::acos(dir.z()/mag(dir));


	if(dir.x() != 0 )
	{ 
		tPhi = Foam::atan(dir.y()/dir.x());
		if(dir.x() < 0 && dir.y() >= 0 ) tPhi = tPhi + pi;
		if(dir.x() < 0 && dir.y() < 0 ) tPhi = tPhi + pi;
		if(dir.x() > 0 && dir.y() < 0 ) tPhi = tPhi + 2*pi;
	}
	else
	{
		if(dir.y() > 0 ) tPhi = pi/2.0;
		if(dir.y() < 0 ) tPhi = 3.0*pi/2.0;;
	}
	
}

void Foam::photoBio::transExteriorSurfaceMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nNbg") << nNbg_ << token::END_STATEMENT << nl;
    os.writeKeyword("nOwn") << nOwn_ << token::END_STATEMENT << nl;
    os.writeKeyword("I0") << I0_ << token::END_STATEMENT << nl;
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamWidthPhi_") << beamWidthPhi_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamDir") << beamDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamNormToSurf") << beamNormToSurf_ << token::END_STATEMENT << nl;
                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{
    makePatchTypeField
    (
        fvPatchScalarField,
        transExteriorSurfaceMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
