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

#include "diffuseEmitterMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "photoBioDOM.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::
diffuseEmitterMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    I0_(0.0), 
    nBands_(1)
{}


Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::
diffuseEmitterMixedFvPatchScalarField
(
    const diffuseEmitterMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    bandDist_(ptf.bandDist_)
{}


Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::
diffuseEmitterMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    I0_(dict.get<scalar>("irradiation")),
    nBands_(dict.get<label>("nBands"))
{
    bandDist_.setSize(nBands_);
    dict.readEntry("bandDist", bandDist_);
}


Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::
diffuseEmitterMixedFvPatchScalarField
(
    const diffuseEmitterMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    bandDist_(ptf.bandDist_)
{}


Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::
diffuseEmitterMixedFvPatchScalarField
(
    const diffuseEmitterMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    I0_(ptf.I0_),
    nBands_(ptf.nBands_),
    bandDist_(ptf.bandDist_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::
updateCoeffs()
{
    // Skip if already updated
    if (this->updated())
    {
        return;
    }
    
    // Change message for parallel communications
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get references to the scalar field and associated models
    scalarField& Iw = *this;
    const photoBioModel& photoBio = db().lookupObject<photoBioModel>("photoBioProperties");
    const photoBioDOM& dom(refCast<const photoBioDOM>(photoBio));
    
    // Ensure model is compatible with boundary condition
    if (dom.nBand() == 0)
    {
        FatalErrorIn
        (
            "Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::updateCoeffs"
        )   << " a non-grey boundary condition is used with a grey"
            << " absorption model" << nl << exit(FatalError);
    }

    // Get the surface normals
    const vectorField n = patch().Sf()/patch().magSf();
    
    // Get the ray id from the field name
    label rayId;
    dom.setRayId(internalField().name(), rayId);

    // Get the band index, ray direction, and solid angle from the ray
    const label iBand = dom.IRay(rayId).iBand();
    const vector& rayDir = dom.IRay(rayId).d();
   
    // Calculate boundary values for all faces
    forAll(Iw, iFace)
    {
        vector surfNorm = -n[iFace];
        scalar cosA = surfNorm & rayDir;

        if (cosA > 0) // Away from the boundary
        {  
            // Boundary value is irradiation per unit solid angle.
            // Since irradiation is specified only away from boundary,
            // value is divided by 2*pi (solid angle for half sphere).
            refValue()[iFace] = I0_*bandDist_[iBand]/(2.0*pi);
            refGrad()[iFace] = 0.0;
            valueFraction()[iFace] = 1.0;
        }
        else // Into the boundary
        {
            refValue()[iFace] = 0.0;
            refGrad()[iFace] = 0.0;
            valueFraction()[iFace] = 0.0;
        }
    }

    // Reset message for parallel communications
    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs(); 
}


void Foam::photoBio::diffuseEmitterMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("irradiation") << I0_ << token::END_STATEMENT << nl;
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("bandDist") << bandDist_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        makePatchTypeField
        (
            fvPatchScalarField,
            diffuseEmitterMixedFvPatchScalarField
        );
    }
}


// ************************************************************************* //
