/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2026 DeGroot Research Group
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "radiationCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

#include "DOM.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::radiationCoupledMixedFvPatchScalarField::
radiationCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    nBands_(1),
    n_(),
    nValidated_(false)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::optical::radiationCoupledMixedFvPatchScalarField::
radiationCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    nBands_(readLabel(dict.lookup("nBands"))),
    n_(),
    nValidated_(false)
{
    dict.lookup("n") >> n_;

    if (n_.size() != nBands_)
    {
        FatalErrorInFunction
            << "n must contain nBands = " << nBands_
            << " entries; got n.size() = " << n_.size()
            << exit(FatalError);
    }

    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "patch type '" << p.type()
            << "' is not a '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        refValue() = scalarField("refValue", dict, p.size());
        refGrad()  = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


Foam::optical::radiationCoupledMixedFvPatchScalarField::
radiationCoupledMixedFvPatchScalarField
(
    const radiationCoupledMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    nBands_(ptf.nBands_),
    n_(ptf.n_),
    nValidated_(ptf.nValidated_)
{}


Foam::optical::radiationCoupledMixedFvPatchScalarField::
radiationCoupledMixedFvPatchScalarField
(
    const radiationCoupledMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    nBands_(ptf.nBands_),
    n_(ptf.n_),
    nValidated_(ptf.nValidated_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::radiationCoupledMixedFvPatchScalarField::validateN_() const
{
    const mappedPatchBase& mpp =
        mappedPatchBase::getMap(patch().patch());
    const label nbrPatchI = mpp.nbrPolyPatch().index();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.nbrMesh());
    const fvPatch& nbrPatch = nbrMesh.boundary()[nbrPatchI];

    const fvPatchScalarField& nbrFieldGeneric =
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            internalField().name()
        );

    // The neighbour MUST be the same BC type; otherwise the user has
    // wired up an asymmetric pair and there's no way to extract its
    // `n` list. Refuse to proceed.
    if (!isA<radiationCoupledMixedFvPatchScalarField>(nbrFieldGeneric))
    {
        FatalErrorInFunction
            << "neighbour patch '" << nbrPatch.name()
            << "' on field '" << internalField().name()
            << "' is type '" << nbrFieldGeneric.type()
            << "', expected 'radiationCoupled'"
            << exit(FatalError);
    }

    const radiationCoupledMixedFvPatchScalarField& nbrField =
        refCast<const radiationCoupledMixedFvPatchScalarField>
        (nbrFieldGeneric);

    if (nbrField.nBands_ != nBands_)
    {
        FatalErrorInFunction
            << "nBands mismatch across radiationCoupled interface: "
            << "this side = " << nBands_
            << ", neighbour = " << nbrField.nBands_
            << exit(FatalError);
    }

    for (label b = 0; b < nBands_; ++b)
    {
        const scalar nOwn = n_[b];
        const scalar nNbg = nbrField.n_[b];
        const scalar relDiff =
            mag(nOwn - nNbg) / max(mag(nOwn), SMALL);
        if (relDiff > nMatchTol)
        {
            FatalErrorInFunction
                << "radiationCoupled requires matching refractive index "
                << "on both sides, but n[" << b << "] differs by "
                << relDiff << " (own = " << nOwn
                << ", neighbour = " << nNbg << ").\n"
                << "    Use 'refractiveCoupled' for genuinely refractive "
                << "interfaces; 'radiationCoupled' is the matched-index "
                << "fast path only."
                << exit(FatalError);
        }
    }

    nValidated_ = true;
}


void Foam::optical::radiationCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!nValidated_)
    {
        validateN_();
    }

    // Tag bump so any processor comms in fromNeighbour() use a
    // distinct message channel from the surrounding solver.
    const int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const DOM& dom = DOM::lookup(db());
    const label rayId = dom.nameToRayId(internalField().name());
    const vector rayDir = dom.IRay(rayId).d();

    // Outward face-normal unit vector (points from owner cell into
    // the boundary). The inward normal is its negative.
    const vectorField nHat = patch().Sf() / patch().magSf();

    // For matched n the cross-region transmission is identity: the
    // same ray d on side A maps to the same ray d on side B with
    // unchanged radiance. So we only need the SAME ray's neighbour-
    // side patch-internal field -- not the full angular sweep that
    // refractiveCoupled does to gather reflected/refracted candidates.
    const mappedPatchBase& mpp =
        mappedPatchBase::getMap(patch().patch());
    const label nbrPatchI = mpp.nbrPolyPatch().index();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.nbrMesh());
    const fvPatch& nbrPatch = nbrMesh.boundary()[nbrPatchI];
    const fvPatchScalarField& nbrField =
        nbrPatch.lookupPatchField<volScalarField, scalar>
        (
            internalField().name()
        );
    const scalarField nbrInternal =
        mpp.fromNeighbour(nbrField.patchInternalField());

    scalarField& Iw = *this;
    forAll(Iw, faceI)
    {
        const vector inward = -nHat[faceI];
        const scalar cosR = rayDir & inward;

        if (cosR > 0.0)
        {
            // Ray flows INTO this domain through this face.
            // Fix the boundary radiance to the neighbour's
            // outgoing-radiance-on-the-same-ray.
            refValue()[faceI] = nbrInternal[faceI];
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
        }
        else
        {
            // Ray flows OUT of this domain. Let the upwind
            // discretisation use the cell-internal value.
            refValue()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 0.0;
        }
    }

    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::optical::radiationCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("n") << n_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{
    makePatchTypeField
    (
        fvPatchScalarField,
        radiationCoupledMixedFvPatchScalarField
    );
} // End namespace optical
} // End namespace Foam

// ************************************************************************* //
