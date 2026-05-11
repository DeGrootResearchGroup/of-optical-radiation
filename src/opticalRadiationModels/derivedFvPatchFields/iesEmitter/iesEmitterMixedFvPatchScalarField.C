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

#include "iesEmitterMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "DOM.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// Floor on cos(d, n_avg) for outgoing rays. Below this, the divergent
// L_d = I/cos relation is replaced by L_d = 0. Picked at 1e-3 so that
// rays within ~3 degrees of grazing are dropped -- well below the
// angular resolution of any DOM grid we run with (nPhi >= 4 gives
// 22.5 deg per cell), so the dropped rays carry negligible flux for
// any IES distribution that is well-behaved at grazing.
static const Foam::scalar IES_MIN_COS = 1e-3;


// Static member definition. See header for rationale.
Foam::HashTable<Foam::scalar>
Foam::optical::iesEmitterMixedFvPatchScalarField::phiTableCache_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::iesEmitterMixedFvPatchScalarField::
iesEmitterMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    iesFile_(),
    nBands_(1),
    fixtureAxis_(0, 0, -1),
    fixtureUp_(1, 0, 0),
    ies_(),
    cacheValid_(false),
    cachedL_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::optical::iesEmitterMixedFvPatchScalarField::
iesEmitterMixedFvPatchScalarField
(
    const iesEmitterMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    iesFile_(ptf.iesFile_),
    nBands_(ptf.nBands_),
    power_(ptf.power_),
    fixtureAxis_(ptf.fixtureAxis_),
    fixtureUp_(ptf.fixtureUp_),
    ies_(),
    cacheValid_(false),
    cachedL_(0.0)
{}


Foam::optical::iesEmitterMixedFvPatchScalarField::
iesEmitterMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    iesFile_(dict.lookup("iesFile")),
    nBands_(readLabel(dict.lookup("nBands"))),
    fixtureAxis_(dict.lookup("fixtureAxis")),
    fixtureUp_(dict.lookup("fixtureUp")),
    ies_(),
    cacheValid_(false),
    cachedL_(0.0)
{
    dict.lookup("power") >> power_;
    if (power_.size() != nBands_)
    {
        FatalErrorInFunction
            << "power must contain nBands = " << nBands_
            << " entries; got " << power_.size()
            << " on patch " << p.name() << " of field " << iF.name()
            << exit(FatalError);
    }
    forAll(power_, b)
    {
        if (power_[b] < 0.0)
        {
            FatalErrorInFunction
                << "power entries must be non-negative; got power[" << b
                << "] = " << power_[b]
                << " on patch " << p.name() << " of field " << iF.name()
                << exit(FatalError);
        }
    }

    // Normalise fixtureAxis_; refuse zero magnitude explicitly so
    // misconfiguration fails loudly rather than producing silent NaN.
    const scalar magAxis = mag(fixtureAxis_);
    if (magAxis < SMALL)
    {
        FatalErrorInFunction
            << "fixtureAxis has zero magnitude on patch " << p.name()
            << " of field " << iF.name()
            << exit(FatalError);
    }
    fixtureAxis_ /= magAxis;

    // Orthogonalise fixtureUp_ against fixtureAxis_, then normalise.
    fixtureUp_ -= (fixtureUp_ & fixtureAxis_)*fixtureAxis_;
    const scalar magUp = mag(fixtureUp_);
    if (magUp < SMALL)
    {
        FatalErrorInFunction
            << "fixtureUp is parallel to fixtureAxis on patch "
            << p.name() << " of field " << iF.name()
            << "; supply a fixtureUp not collinear with fixtureAxis"
            << exit(FatalError);
    }
    fixtureUp_ /= magUp;
}


Foam::optical::iesEmitterMixedFvPatchScalarField::
iesEmitterMixedFvPatchScalarField
(
    const iesEmitterMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    iesFile_(ptf.iesFile_),
    nBands_(ptf.nBands_),
    power_(ptf.power_),
    fixtureAxis_(ptf.fixtureAxis_),
    fixtureUp_(ptf.fixtureUp_),
    ies_(),
    cacheValid_(false),
    cachedL_(0.0)
{}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::scalar
Foam::optical::iesEmitterMixedFvPatchScalarField::gammaDegFromDir_
(
    const vector& d
) const
{
    const scalar cosG = max(-1.0, min(1.0, d & fixtureAxis_));
    return Foam::acos(cosG)*180.0/pi;
}


Foam::scalar
Foam::optical::iesEmitterMixedFvPatchScalarField::hDegFromDir_
(
    const vector& d
) const
{
    // Project d into the plane perpendicular to fixtureAxis_.
    const vector dPerp = d - (d & fixtureAxis_)*fixtureAxis_;
    const scalar magPerp = mag(dPerp);
    if (magPerp < SMALL)
    {
        // d is along ±fixtureAxis_; horizontal angle is undefined,
        // pick 0 (the IES table value at the pole is independent of
        // h in well-formed files).
        return 0.0;
    }
    // Right-handed local frame on the perpendicular plane:
    //   e1 = fixtureUp_   (h = 0)
    //   e2 = fixtureAxis_ x fixtureUp_   (h = 90)
    const vector e2 = fixtureAxis_ ^ fixtureUp_;
    const scalar cosH = (dPerp & fixtureUp_)/magPerp;
    const scalar sinH = (dPerp & e2)/magPerp;
    scalar h = Foam::atan2(sinH, cosH)*180.0/pi;
    if (h < 0.0) h += 360.0;
    return h;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::iesEmitterMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    scalarField& Iw = *this;
    const DOM& dom = DOM::lookup(db());

    if (dom.nBand() != nBands_)
    {
        FatalErrorInFunction
            << "nBands in BC (" << nBands_ << ") does not match the model's"
            << " nBand (" << dom.nBand() << ")"
            << exit(FatalError);
    }

    // This BC instance is bound to one specific I_<band>_<angle>
    // field. Find that ray's id, band and direction.
    const label rayId = dom.nameToRayId(internalField().name());
    const label iBand = dom.IRay(rayId).iBand();
    const vector& rayDir = dom.IRay(rayId).d();

    if (!cacheValid_)
    {
        // Lazy IES load: deferring to first updateCoeffs avoids
        // touching the disk during BC construction (which can run in
        // surprising contexts -- mapping, decomposition, etc.).
        if (!ies_.valid())
        {
            ies_.reset(new iesPhotometry(iesFile_));
        }

        // Global patch geometry. magSf and Sf are local; reduce across
        // processors so the renormalisation is patch-global rather
        // than per-rank.
        scalar globalArea = sum(patch().magSf());
        vector globalSf = sum(patch().Sf());
        reduce(globalArea, sumOp<scalar>());
        reduce(globalSf, sumOp<vector>());

        if (globalArea < SMALL)
        {
            FatalErrorInFunction
                << "Patch " << patch().name() << " has zero area"
                << exit(FatalError);
        }

        // Patch-averaged inward normal (into the domain). Sf points
        // out of the domain by OpenFOAM convention, so negate.
        const vector nAvgInDomain = -globalSf/globalArea;

        // Phi_table: sum over rays in this band that go INTO the
        // domain through the patch (cos > floor) of I_table(d)*Omega.
        // No cos weighting: this normalisation makes
        //   sum_d L_d * Omega_d * A_proj(d) == P
        // with A_proj(d) = A_patch * (d.n_avg)+, exactly when L_d is
        // formed as I_table/cos (see header derivation).
        //
        // Phi_table is identical across all nAngle ray-BC instances
        // on a given (patch, band), so consult the cross-instance
        // cache first. Key is patch name + band index. The cache is
        // populated lazily on the first BC instance to reach this
        // point for a given (patch, band); subsequent instances reuse.
        const label nAngle = dom.nAngle();
        const word cacheKey = patch().name() + "_b" + Foam::name(iBand);
        scalar phiTable;
        if (phiTableCache_.found(cacheKey))
        {
            phiTable = phiTableCache_[cacheKey];
        }
        else
        {
            phiTable = 0.0;
            for (label iAngle = 0; iAngle < nAngle; ++iAngle)
            {
                const ray& rj = dom.IRay(iBand*nAngle + iAngle);
                const scalar cosIn = rj.d() & nAvgInDomain;
                if (cosIn > IES_MIN_COS)
                {
                    const scalar Ij = ies_->interpolate
                    (
                        gammaDegFromDir_(rj.d()),
                        hDegFromDir_(rj.d())
                    );
                    phiTable += Ij*rj.omega();
                }
            }
            phiTableCache_.insert(cacheKey, phiTable);
        }

        // This ray's contribution.
        const scalar cosThisRay = rayDir & nAvgInDomain;
        if (cosThisRay > IES_MIN_COS && phiTable > VSMALL)
        {
            const scalar I = ies_->interpolate
            (
                gammaDegFromDir_(rayDir),
                hDegFromDir_(rayDir)
            );
            cachedL_ =
                (power_[iBand]/(globalArea*phiTable))
               *I/cosThisRay;
        }
        else
        {
            cachedL_ = 0.0;
        }

        cacheValid_ = true;
    }

    const vectorField n = patch().Sf()/patch().magSf();

    forAll(Iw, iFace)
    {
        // Inward-pointing surface normal (into the domain).
        const vector surfNorm = -n[iFace];
        const scalar cosFace = surfNorm & rayDir;

        if (cosFace > 0.0)
        {
            // Ray going INTO the domain through this face; pin the
            // value (fixedValue behaviour from the mixed BC).
            valueFraction()[iFace] = 1.0;
            refValue()[iFace] = cachedL_;
            refGrad()[iFace] = 0.0;
        }
        else
        {
            // Ray going OUT of the domain through this face;
            // zeroGradient so the cell radiance flows through.
            valueFraction()[iFace] = 0.0;
            refValue()[iFace] = 0.0;
            refGrad()[iFace] = 0.0;
        }
    }

    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::optical::iesEmitterMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("iesFile") << iesFile_ << token::END_STATEMENT << nl;
    os.writeKeyword("power") << power_ << token::END_STATEMENT << nl;
    os.writeKeyword("fixtureAxis") << fixtureAxis_
        << token::END_STATEMENT << nl;
    os.writeKeyword("fixtureUp") << fixtureUp_
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
            iesEmitterMixedFvPatchScalarField
        );
    }
}


// ************************************************************************* //
