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

#include "surroundingMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "DOM.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    surroundOptical_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::optical::surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const surroundingMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    surroundOptical_(ptf.surroundOptical_)
{}


Foam::optical::surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    surroundOptical_(readScalar(dict.lookup("surroundOptical")))
{


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


Foam::optical::surroundingMixedFvPatchScalarField::
surroundingMixedFvPatchScalarField
(
    const surroundingMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    surroundOptical_(ptf.surroundOptical_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::surroundingMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
 // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

      
    scalarField& Iw = *this;
    
    const radiationModel& opticalRadiation = db().lookupObject<radiationModel>("opticalRadiationProperties");

    const DOM& dom(refCast<const DOM>(opticalRadiation));
 
    
    label rayId = -1;
    dom.setRayId(internalField().name(), rayId);

    const  vectorField n = patch().Sf()/patch().magSf();
    
    const  vector& bdRayDir = dom.IRay(rayId).d();
    const  scalar bdOmega = dom.IRay(rayId).omega();
    const  label nBand = dom.nBand();  
	
    forAll(Iw, faceI)
    {
	
	   if((-n[faceI] & bdRayDir) > 0.0 )    // direction out of the wall   
       {	
            refValue()[faceI] = surroundOptical_/nBand/4/pi ; ; 
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


void Foam::optical::surroundingMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("surroundOptical  ")  << surroundOptical_ << token::END_STATEMENT << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{
    makePatchTypeField
    (
        fvPatchScalarField,
        surroundingMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
