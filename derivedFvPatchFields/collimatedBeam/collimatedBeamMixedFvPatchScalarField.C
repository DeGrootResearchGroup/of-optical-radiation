/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

#include "collimatedBeamMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "DOM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::collimatedBeamMixedFvPatchScalarField::
collimatedBeamMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    beamDirection_(vector::zero),
    nBands_(1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::optical::collimatedBeamMixedFvPatchScalarField::
collimatedBeamMixedFvPatchScalarField
(
    const collimatedBeamMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    beamDirection_(ptf.beamDirection_),
    nBands_(ptf.nBands_),
    beamRadiance_(ptf.beamRadiance_)
{}


Foam::optical::collimatedBeamMixedFvPatchScalarField::
collimatedBeamMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    beamDirection_(dict.lookup("beamDirection")),
    nBands_(readLabel(dict.lookup("nBands")))
{
    // Normalise the beam direction; refuse a zero vector explicitly so
    // misconfiguration fails loudly rather than silently producing a NaN
    // direction lookup.
    const scalar magBeam = mag(beamDirection_);
    if (magBeam < SMALL)
    {
        FatalErrorInFunction
            << "beamDirection has zero magnitude on patch " << p.name()
            << " of field " << iF.name()
            << exit(FatalError);
    }
    beamDirection_ /= magBeam;

    dict.lookup("beamRadiance") >> beamRadiance_;
    if (beamRadiance_.size() != nBands_)
    {
        FatalErrorInFunction
            << "beamRadiance must contain nBands = " << nBands_
            << " entries; got " << beamRadiance_.size()
            << exit(FatalError);
    }
}


Foam::optical::collimatedBeamMixedFvPatchScalarField::
collimatedBeamMixedFvPatchScalarField
(
    const collimatedBeamMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    beamDirection_(ptf.beamDirection_),
    nBands_(ptf.nBands_),
    beamRadiance_(ptf.beamRadiance_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::collimatedBeamMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    scalarField& Iw = *this;
    const radiationModel& opticalRadiation =
        db().lookupObject<radiationModel>("opticalRadiationProperties");
    const DOM& dom(refCast<const DOM>(opticalRadiation));

    if (dom.nBand() != nBands_)
    {
        FatalErrorInFunction
            << "nBands in BC (" << nBands_ << ") does not match the model's"
            << " nBand (" << dom.nBand() << ")"
            << exit(FatalError);
    }

    // This BC is instantiated for one specific I_<band>_<angle> field.
    // Find that ray's id and band.
    const label rayId = dom.nameToRayId(internalField().name());
    const label iBand = dom.IRay(rayId).iBand();
    const vector& rayDir = dom.IRay(rayId).d();

    // Determine which discrete ray bin contains the beam direction in this
    // band. Only that ray sees a non-zero refValue from the beam.
    const label beamRayId = dom.dirToRayId(beamDirection_, iBand);
    const bool isBeamRay = (rayId == beamRayId);

    const vectorField n = patch().Sf()/patch().magSf();

    forAll(Iw, iFace)
    {
        // Inward-pointing surface normal (into the domain).
        const vector surfNorm = -n[iFace];
        const scalar cosRay  = surfNorm & rayDir;

        if (cosRay > 0.0)
        {
            // This ray is going INTO the domain through this face --
            // pin its value (fixedValue behaviour from the mixed BC).
            valueFraction()[iFace] = 1.0;
            refGrad()[iFace] = 0.0;

            if (isBeamRay)
            {
                // Beam only illuminates faces whose inward normal has a
                // positive component along the beam direction.
                const scalar cosBeam = surfNorm & beamDirection_;
                refValue()[iFace] =
                    (cosBeam > 0.0) ? beamRadiance_[iBand] : 0.0;
            }
            else
            {
                refValue()[iFace] = 0.0;
            }
        }
        else
        {
            // Ray going OUT of the domain -- zeroGradient so the cell
            // value flows through.
            valueFraction()[iFace] = 0.0;
            refValue()[iFace] = 0.0;
            refGrad()[iFace] = 0.0;
        }
    }

    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::optical::collimatedBeamMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("beamDirection") << beamDirection_
        << token::END_STATEMENT << nl;
    os.writeKeyword("beamRadiance") << beamRadiance_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        makePatchTypeField
        (
            fvPatchScalarField,
            collimatedBeamMixedFvPatchScalarField
        );
    }
}


// ************************************************************************* //
