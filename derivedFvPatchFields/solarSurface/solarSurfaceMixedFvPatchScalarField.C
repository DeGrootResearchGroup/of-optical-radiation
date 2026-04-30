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

#include "solarSurfaceMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "DOM.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::solarSurfaceMixedFvPatchScalarField::
solarSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(0.0),
    nOwn_(0.0), 
    diffuseCoeff_(0.0),

    directBeam_(0.0), 
    skyDiffuse_(0.0), 
    groundReflected_(0.0), 
    zenithAngle_(0.0),
    azimuthAngle_(0.0),
 //   visualAngle_(0.0),
    
    nBands_(1)
 //   opticalRadiationBandDist_(null),

{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::optical::solarSurfaceMixedFvPatchScalarField::
solarSurfaceMixedFvPatchScalarField
(
    const solarSurfaceMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    diffuseCoeff_(ptf.diffuseCoeff_), 
    directBeam_(ptf.directBeam_), 
    skyDiffuse_(ptf.skyDiffuse_), 
    groundReflected_(ptf.groundReflected_), 
    zenithAngle_(ptf.zenithAngle_), 
    azimuthAngle_(ptf.azimuthAngle_), 
 //   visualAngle_(ptf.visualAngle_), 

    nBands_(ptf.nBands_)
 //   opticalRadiationBandDist_(ptf.opticalRadiationBandDist_)
{}


Foam::optical::solarSurfaceMixedFvPatchScalarField::
solarSurfaceMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nNbg_(readScalar(dict.lookup("nNbg"))), //scalarField("n1", dict),
    nOwn_(readScalar(dict.lookup("nOwn"))), // scalarField("n2", dict);
    diffuseCoeff_(readScalar(dict.lookup("diffuseCoeff"))),
    directBeam_(readScalar(dict.lookup("directBeam"))),
    skyDiffuse_(readScalar(dict.lookup("skyDiffuse"))),
    groundReflected_(readScalar(dict.lookup("groundReflected"))),
    zenithAngle_(readScalar(dict.lookup("zenithAngle"))),
    azimuthAngle_(readScalar(dict.lookup("azimuthAngle"))),
//    visualAngle_(readScalar(dict.lookup("visualAngle"))),
    nBands_(readLabel(dict.lookup("nBand"))) // scalarField("n2", dict);
{
	
   opticalRadiationBandDist_.setSize(nBands_);   
   dict.lookup("opticalRadiationBandDist") >> opticalRadiationBandDist_;
	  
    if (dict.found("refValue"))
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



Foam::optical::solarSurfaceMixedFvPatchScalarField::
solarSurfaceMixedFvPatchScalarField
(
    const solarSurfaceMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_),
    diffuseCoeff_(ptf.diffuseCoeff_), 
    directBeam_(ptf.directBeam_), 
    skyDiffuse_(ptf.skyDiffuse_), 
    groundReflected_(ptf.groundReflected_), 
    zenithAngle_(ptf.zenithAngle_), 
    azimuthAngle_(ptf.azimuthAngle_), 
  //  visualAngle_(ptf.visualAngle_), 
    nBands_(ptf.nBands_)
{}






// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::solarSurfaceMixedFvPatchScalarField::
updateCoeffs()
{
	if (this->updated())
    {
        return;
    }
    
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

     	
    scalarField& Iw = *this;
   
    const radiationModel& opticalRadiation = db().lookupObject<radiationModel>("opticalRadiationProperties");

    const DOM& dom(refCast<const DOM>(opticalRadiation));
 
    if (dom.nBand() == 0)
    {
        FatalErrorIn
        (
            "Foam::optical::"
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
 //   const  vector& bdRayDAve = dom.IRay(rayId).dAve();
    const  scalar bdOmega = dom.IRay(rayId).omega();
    const  scalar bdRayPhi  = dom.IRay(rayId).phi();
    const  scalar bdRayTheta = dom.IRay(rayId).theta();

    const  label nTheta = dom.nTheta();    
    const  label nPhi = dom.nPhi();    
    const  scalar deltaPhi   =  pi /(2.0*nPhi);
    const  scalar deltaTheta =  pi  /nTheta;
    
    label  npPhi  = dom.nPixelPhi();    
    label  npTheta = dom.nPixelTheta();  

	scalar  nRatio = nOwn_/  nNbg_;  
	
    const scalar sinTheta = Foam::sin((180-zenithAngle_)/180*pi);
    const scalar cosTheta = Foam::cos((180-zenithAngle_)/180*pi);
    const scalar sinPhi = Foam::sin((azimuthAngle_-180)/180*pi);
    const scalar cosPhi = Foam::cos((azimuthAngle_-180)/180*pi);
   
    const vector solarDir = vector(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
    label  solarRay = -1;
    dom.dirToRayId(solarDir, 0, solarRay);	// this is important 
    
	scalar specularReflection = 0.0;
	scalar diffusiveReflection =0.0;
	scalar specularRefraction = 0.0;
	scalar diffusiveRefraction =0.0;
	scalar R = 0.0;
	

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
       specularRefraction = 0.0;
       diffusiveRefraction = 0.0;
       
	   vector surfNorm = -n[faceI] ;
	   scalar cosB = surfNorm & bdRayDir; 
	              
		if ( cosB > 0.0)   
		{
	
			for (label jAngle = 0; jAngle < nAngle; jAngle++)
			{         
				label sweepRayID = jAngle + (iBand)*nAngle;
				vector sweepDir = dom.IRay(sweepRayID).d();
				vector sweepdAve = dom.IRay(sweepRayID).dAve();   
				
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
					    
						scalar r1 = (nNbg_*cosB - nOwn_*cosA ) / (nNbg_*cosB + nOwn_*cosA );
						scalar r2 = (nNbg_*cosA - nOwn_*cosB ) / (nNbg_*cosA + nOwn_*cosB );
            
					    R = 0.5*(r1*r1 + r2*r2); 
						
          			    label refracIncidentRay = -1;  
          		        dom.dirToRayId(refracIncidentDir, 0, refracIncidentRay);	// this is important 
						
						if(refracIncidentDir[2] < 0)
						{
							diffusiveRefraction = diffusiveRefraction + skyDiffuse_/2/pi*(1-R)*mag(surfNorm & sweepdAve) ; ;  
						}
						else
						{
							diffusiveRefraction = diffusiveRefraction + groundReflected_/2/pi*(1-R)*mag(surfNorm & sweepdAve) ;  
						}
						
						if(refracIncidentRay == solarRay)
						{
							diffusiveRefraction = diffusiveRefraction +  directBeam_*(1-R)*mag(surfNorm & sweepdAve) ;  
						}
					 }
					 else
					 {
						 R = 1;
					 }
					 
						diffusiveReflection = diffusiveReflection + reflecFace[faceI]*R*mag(surfNorm & sweepdAve) ;  
                        
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

                    
                    if( cosB*cosB > 1-1/(nRatio*nRatio) )    
					{ 
                    
						vector refracIncidentDir = (pixelDir - cosB * surfNorm )*nRatio 
						+ Foam::sqrt(1 -(nRatio*nRatio)*(1- cosB*cosB))*surfNorm ;
					
						scalar cosA = surfNorm& pixelDir; // / mag(d);                //    /mag(n[faceI])
          
						scalar r1 = (nNbg_*cosB - nOwn_*cosA ) / (nNbg_*cosB + nOwn_*cosA );
						scalar r2 = (nNbg_*cosA - nOwn_*cosB ) / (nNbg_*cosA + nOwn_*cosB );
            
						scalar R =  0.5*(r1*r1 + r2*r2);
   
          			    label refracIncidentRay = -1;  
          		        dom.dirToRayId(refracIncidentDir, 0, refracIncidentRay);	// this is important 
						
						if(refracIncidentDir[2] < 0)
						{
							specularRefraction = specularRefraction + skyDiffuse_/pi/2*(1-R)*pixelOmega ; ;  
						}
						else
						{
							specularRefraction = specularRefraction + groundReflected_/2/pi*(1-R)*pixelOmega ;  
						}
						
						if(refracIncidentRay == solarRay)
						{
							specularRefraction = specularRefraction +  directBeam_*(1-R)*pixelOmega ;  
						}
					}
					else
					{
						R = 1.0 ;
					}
					
					specularReflection = specularReflection + reflecFace[faceI]*R*pixelOmega;   
				}
			}
		}
           
                          
    //       refValue()[faceI] = (diffuseCoeff_*(diffusiveRefraction+ diffusiveReflection )/pi/2 + (1.0 - diffuseCoeff_)*(specularReflection + specularRefraction)/bdOmega)/bdOmega;  //careful here, first omega convert the opticalRadiation transmision control angle, second one convert the value to angle intensity 
            refValue()[faceI] = (diffuseCoeff_*diffusiveRefraction/pi/2 + (1.0 - diffuseCoeff_)*specularRefraction/bdOmega)/bdOmega;     
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

void Foam::optical::solarSurfaceMixedFvPatchScalarField::dirToAngle
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


void Foam::optical::solarSurfaceMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nNbg") << nNbg_ << token::END_STATEMENT << nl;
    os.writeKeyword("nOwn") << nOwn_ << token::END_STATEMENT << nl;
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("directBeam") << directBeam_ << token::END_STATEMENT << nl;
    os.writeKeyword("skyDiffuse") << skyDiffuse_ << token::END_STATEMENT << nl;
    os.writeKeyword("groundReflected") << groundReflected_ << token::END_STATEMENT << nl;
    os.writeKeyword("zenithAngle") << zenithAngle_ << token::END_STATEMENT << nl;
    os.writeKeyword("azimuthAngle") << azimuthAngle_ << token::END_STATEMENT << nl;
   // os.writeKeyword("visualAngle") << visualAngle_ << token::END_STATEMENT << nl;
                
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{
    makePatchTypeField
    (
        fvPatchScalarField,
        solarSurfaceMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
