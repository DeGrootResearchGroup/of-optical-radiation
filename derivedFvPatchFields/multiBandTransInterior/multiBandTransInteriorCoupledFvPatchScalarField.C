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

#include "multiBandTransInteriorCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

#include "DOM.H"
#include "constants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::multiBandTransInteriorCoupledFvPatchScalarField::
multiBandTransInteriorCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    diffuseFraction_(0.0),
    nBands_(1)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


Foam::optical::multiBandTransInteriorCoupledFvPatchScalarField::
multiBandTransInteriorCoupledFvPatchScalarField
(
    const multiBandTransInteriorCoupledFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    diffuseFraction_(ptf.diffuseFraction_),
    nBands_(ptf.nBands_),
    nNbg_(ptf.nNbg_),
    nOwn_(ptf.nOwn_)
{}


Foam::optical::multiBandTransInteriorCoupledFvPatchScalarField::
multiBandTransInteriorCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    mixedFvPatchScalarField(p, iF),
    diffuseFraction_(readScalar(dict.lookup("diffuseFraction"))),
    nBands_(readLabel(dict.lookup("nBands")))
{
    // Read the refractive indices for each band; the stream extractor
    // sizes the lists, so we just verify the size matches nBands_.
    dict.lookup("nNbg") >> nNbg_;
    dict.lookup("nOwn") >> nOwn_;

    if (nNbg_.size() != nBands_ || nOwn_.size() != nBands_)
    {
        FatalErrorInFunction
            << "nNbg and nOwn must each contain nBands = " << nBands_
            << " entries; got nNbg.size() = " << nNbg_.size()
            << ", nOwn.size() = " << nOwn_.size()
            << exit(FatalError);
    }

    // Check that the patch is of the correct type
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
            (
                "multiBandTransInteriorCoupledFvPatchScalarField::"
                "multiBandTransInteriorCoupledFvPatchScalarField\n"
                "(\n"
                "    const fvPatch& p,\n"
                "    const DimensionedField<scalar, volMesh>& iF,\n"
                "    const dictionary& dict\n"
                ")\n"
            )   << "\n    patch type '" << p.type()
                << "' not type '" << mappedPatchBase::typeName << "'"
                << "\n    for patch " << p.name()
                << " of field " << internalField().name()
                << " in file " << internalField().objectPath()
                << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad()  = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {

        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;

    }
}


Foam::optical::multiBandTransInteriorCoupledFvPatchScalarField::
multiBandTransInteriorCoupledFvPatchScalarField
(
    const multiBandTransInteriorCoupledFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    mixedFvPatchScalarField(wtcsf, iF),
    diffuseFraction_(wtcsf.diffuseFraction_),
    nBands_(wtcsf.nBands_),
    nNbg_(wtcsf.nNbg_),
    nOwn_(wtcsf.nOwn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::multiBandTransInteriorCoupledFvPatchScalarField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the cell and boundary patch values
    scalarField Ic(patchInternalField());
    scalarField& Iw = *this;

    // Get the opticalRadiation model and the associated DOM object
    const radiationModel& opticalRadiation = db().lookupObject<radiationModel>("opticalRadiationProperties");
    const DOM& dom(refCast<const DOM>(opticalRadiation));

    // Get the ray ID for this ray
    label rayId = dom.nameToRayId(internalField().name());

    // Get the surface normals
    const label patchI = patch().index();
    vectorField n = patch().Sf()/patch().magSf();

    // Get the angular discretization information
    const label iBand = dom.IRay(rayId).iBand();
    const vector rayDir = dom.IRay(rayId).d();
    const scalar deltaTheta = dom.deltaTheta();
    const scalar deltaPhi = dom.deltaPhi();
    const label nAngle = dom.nAngle();
    const scalar rayPhi = dom.IRay(rayId).phi();
    const scalar rayTheta = dom.IRay(rayId).theta();

    // Get the pixel discretization information
    label npPhi = dom.nPixelPhi();
    label npTheta = dom.nPixelTheta();

    // Correct the pixel discretization information if case is 2D
    if (internalField().mesh().nSolutionD() == 2)
    {
        npTheta = 1;
    }

    // Create a list of scalarField objects to store the other rays
    PtrList<scalarField> NbrRaySet(nAngle*nBands_);

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label nbrPatchI = mpp.nbrPolyPatch().index();
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.nbrMesh());
    const fvPatch& nbrPatch = nbrMesh.boundary()[nbrPatchI];

    // Loop through all rays and store the internal cell values for the patch
    for (label jBand = 0; jBand < nBands_; jBand++)
    {
        for (label jAngle = 0; jAngle < nAngle; jAngle++)
        {
            label rayI = jAngle + jBand*nAngle;

            const multiBandTransInteriorCoupledFvPatchScalarField&
                nbrField = refCast
                <const multiBandTransInteriorCoupledFvPatchScalarField>
                (nbrPatch.lookupPatchField<volScalarField, scalar>(dom.IRay(rayI).I().name()));

            NbrRaySet.set(rayI, mpp.fromNeighbour(nbrField.patchInternalField()));
        }
    }

    // Store refractive indices for each side, as well as their ratio
    const scalar nA = nNbg_[iBand];
    const scalar nB = nOwn_[iBand];
    const scalar nRatio = nB/nA;

    // Loop through all faces on the boundary patch
    forAll(Iw, faceI)
    {
        // Initialize values for specular and diffuse radiation for this face
        scalar specular = 0.0;
        scalar diffuse = 0.0;

        // Initialize sums of integral of direction vector over solid angle
        scalar magIntDirOmega = 0.0;
        scalar posIntDirOmega = 0.0;

        // Get the surface normal for this face (directed into domain)
        vector surfNorm = -n[faceI];

        // Calculate cosine of angle between this ray and surface normal
        scalar cosR = rayDir & surfNorm;

        // Check for control angle overhang
        bool overhang = false;
        for (label i = 1; i <= npTheta; i++)
        {
            scalar pxRayTheta = rayTheta - 0.5*deltaTheta + (i - 0.5)*deltaTheta/npTheta;
            for (label j = 1; j <= npPhi; j++)
            {
                scalar pxRayPhi = rayPhi - 0.5*deltaPhi + (j - 0.5)*deltaPhi/npPhi;
                vector pixelDir = dom.anglesToDir(pxRayTheta, pxRayPhi);
                if (cosR*(pixelDir & surfNorm) < 0.0)
                {
                    overhang = true;
                }
            }
        }

        // Case 1: Ray is coming out of wall, or there is overhang -> Calculate intensity at the wall
        if (cosR > 0.0 || overhang)
        {
            // Note: diffuse radiation code has not been checked at all; use at own risk!
            if (diffuseFraction_ > 0)
            {
                for (label jAngle = 0; jAngle < nAngle; jAngle++)
                {
                    label sweepRayId = jAngle + iBand*nAngle;

                    vector sweepDir = dom.IRay(sweepRayId).d();
                    vector sweepdAve = dom.IRay(sweepRayId).dAve();

                    scalar cosB = surfNorm & sweepDir;

                    if (cosB > 0.0) // direction out of the wall
                    {
                        vector reflecIncidentDir = sweepDir - 2*cosB*surfNorm;
                        label reflecIncidentRay = dom.dirToRayId(reflecIncidentDir, iBand);
                        const scalarField& reflecFace = dom.IRay(reflecIncidentRay).I().boundaryField()[patchI];

                        if (cosB*cosB > 1 - 1/(nRatio*nRatio))
                        {
                            vector refracIncidentDir = (sweepDir - cosB*surfNorm)*nRatio
                                + Foam::sqrt(1 -(nRatio*nRatio)*(1- cosB*cosB))*surfNorm;

                            scalar cosA = surfNorm & refracIncidentDir;

                            scalar r1 = (nA*cosB - nB*cosA)/(nA*cosB + nB*cosA);
                            scalar r2 = (nA*cosA - nB*cosB)/(nA*cosA + nB*cosB);
                            scalar R = 0.5*(r1*r1 + r2*r2);

                            label refracIncidentRay = dom.dirToRayId(refracIncidentDir, iBand);

                            // Transmitted radiance picks up an (nB/nA)^2 factor
                            // from solid-angle compression across the interface
                            // (radiance invariant is I/n^2, not I).
                            diffuse += (NbrRaySet[refracIncidentRay][faceI]*nRatio*nRatio*(1-R) + reflecFace[faceI]*R) *mag(surfNorm & sweepdAve);

                        }
                        else
                        {
                            diffuse = diffuse + reflecFace[faceI]*mag(surfNorm & sweepdAve);
                        }
                    }
                }
            }

            // Create a list to store candidate rays for reflection and refraction
            DynamicList<label> candidateReflectedRays;
            DynamicList<label> candidateRefractedRays;

            // Determine the rays contributing to outgoing intensity, looping though each theta pixel
            for (label i = 1; i <= npTheta; i++)
            {
                // Calculate the theta angle for this pixel
                scalar pxRayTheta = rayTheta - 0.5*deltaTheta + (i - 0.5)*deltaTheta/npTheta;

                // Loop through all phi pixels
                for (label j = 1; j <= npPhi; j++)
                {
                    // Calculate the phi angle for this pixel
                    scalar pxRayPhi = rayPhi - 0.5*deltaPhi + (j - 0.5)*deltaPhi/npPhi;

                    // Create a vector in the direction of this pixel
                    vector pixelDir = dom.anglesToDir(pxRayTheta, pxRayPhi);

                    // Get the cosine of the angle between the pixel direction and the surface normal
                    scalar cosB = pixelDir & surfNorm;

                    // Calculate integral of direction vector wrt solid angle over pixel
                    vector intDirOmega = dom.intDirOmega(pxRayTheta, pxRayPhi, deltaTheta/npTheta, deltaPhi/npPhi);

                    // Pixel ray is coming out of the surface (store incoming rays)
                    if (cosB > 0.0)
                    {
                        // Calculate direction of reflected incident radiation for outgoing pixel ray
                        vector reflectIncidentDir = pixelDir - 2.0*cosB*surfNorm;

                        // Get the ID of the ray in the direction of the reflected incident ray
                        label reflectIncidentRay = dom.dirToRayId(reflectIncidentDir, iBand);

                        // Add this ray to the list of candidate rays if it is not already added
                        if (std::find(candidateReflectedRays.begin(), candidateReflectedRays.end(), reflectIncidentRay) == candidateReflectedRays.end())
                        {
                            candidateReflectedRays.append(reflectIncidentRay);
                        }

                        // Check that outgoing pixel ray is within the critical angle for refracted rays coming from side A
                        if ((1 - (nRatio*nRatio)*(1 - cosB*cosB)) > 0.0)
                        {
                            // Calculate direction of refracted radiation
                            vector refractIncidentDir = (pixelDir - cosB*surfNorm)*nRatio + Foam::sqrt(1 - (nRatio*nRatio)*(1 - cosB*cosB))*surfNorm;

                            // Get the ID of the ray in the direction of the refracted incident ray
                            label refractIncidentRay = dom.dirToRayId(refractIncidentDir, iBand);

                            // Add this ray to the list of candidate rays if it is not already added
                            if (std::find(candidateRefractedRays.begin(), candidateRefractedRays.end(), refractIncidentRay) == candidateRefractedRays.end())
                            {
                                candidateRefractedRays.append(refractIncidentRay);
                            }
                        }

                        // Accumulate integrals of direction vector over solid angle
                        magIntDirOmega += intDirOmega & surfNorm;
                        posIntDirOmega += intDirOmega & surfNorm;
                    }
                    // Pixel ray is going into the surface (accumulate incident energy)
                    else
                    {
                        // Accumulate integrals of direction vector over solid angle
                        magIntDirOmega -= intDirOmega & surfNorm;
                    }
                }
            }

            // Calculate specularly reflected energy, looping through all candidate rays
            forAll(candidateReflectedRays, rayIdx)
            {
                // Calculate the ray angles
                label iRay = candidateReflectedRays[rayIdx];
                scalar iRayTheta = dom.IRay(iRay).theta();
                scalar iRayPhi = dom.IRay(iRay).phi();

                // Loop though all theta pixels
                for (label i = 1; i <= npTheta; i++)
                {
                    // Calculate the theta angle for this pixel
                    scalar pxRayTheta = iRayTheta - 0.5*deltaTheta + (i - 0.5)*deltaTheta/npTheta;

                    // Loop through all phi pixels
                    for (label j = 1; j <= npPhi; j++)
                    {
                        // Calculate the phi angle for this pixel
                        scalar pxRayPhi = iRayPhi - 0.5*deltaPhi + (j - 0.5)*deltaPhi/npPhi;

                        // Create a vector in the direction of this pixel
                        vector pixelDir = dom.anglesToDir(pxRayTheta, pxRayPhi);

                        // Get the cosine of the angle between the pixel direction and the surface normal
                        scalar cosP = - pixelDir & surfNorm;

                        // Create a vector in the direction of the reflected ray
                        vector reflectDir = pixelDir + 2.0*cosP*surfNorm;

                        // Cosine of the outgoing reflected ray with the
                        // surface normal. reflectDir is unit-length by
                        // construction, so the dot product gives the exact
                        // continuous cosine -- Fresnel R is a physical
                        // relationship between continuous angles, not the
                        // discrete pixel grid.
                        scalar cosB = reflectDir & surfNorm;

                        // Continue only if reflected ray is coming out of the surface
                        if (cosB > 0.0)
                        {
                            // Initialize reflectivity for total internal reflection
                            scalar R = 1.0;

                            // Check if the the outgoing ray results from refraction
                            if ((1 - (nRatio*nRatio)*(1 - cosB*cosB)) > 0.0)
                            {
                                // Calculate direction of refracted radiation
                                vector refractIncidentDir = (reflectDir - cosB*surfNorm)*nRatio + Foam::sqrt(1 - (nRatio*nRatio)*(1 - cosB*cosB))*surfNorm;

                                // Cosine of the refracted incident ray with the
                                // surface normal -- exact continuous value.
                                scalar cosA = refractIncidentDir & surfNorm;

                                // Calculate interface reflectivity
                                scalar r1 = (nA*cosB - nB*cosA)/(nA*cosB + nB*cosA);
                                scalar r2 = (nA*cosA - nB*cosB)/(nA*cosA + nB*cosB);
                                R = 0.5*(r1*r1 + r2*r2);
                            }

                            // Calculate integral of direction vector wrt solid angle over pixel
                            vector intDirOmega = dom.intDirOmega(pxRayTheta, pxRayPhi, deltaTheta/npTheta, deltaPhi/npPhi);

                            // Get the ID of the ray in the direction of the outgoing reflected ray
                            label reflectRay = dom.dirToRayId(reflectDir, iBand);

                            // Accumulate the specular energy only if the outgoing ray matches this ray
                            if (reflectRay == rayId)
                            {
                                specular += R*dom.IRay(iRay).I().boundaryField()[patchI][faceI]*(-intDirOmega & surfNorm);
                            }
                        }
                    }
                }
            }

            // Calculate specularly refracted energy, looping through all candidate rays
            forAll(candidateRefractedRays, rayIdx)
            {
                // Calculate the ray angles
                label iRay = candidateRefractedRays[rayIdx];
                scalar iRayTheta = dom.IRay(iRay).theta();
                scalar iRayPhi = dom.IRay(iRay).phi();

                // Loop though all theta pixels
                for (label i = 1; i <= npTheta; i++)
                {
                    // Calculate the theta angle for this pixel
                    scalar pxRayTheta = iRayTheta - 0.5*deltaTheta + (i - 0.5)*deltaTheta/npTheta;

                    // Loop through all phi pixels
                    for (label j = 1; j <= npPhi; j++)
                    {
                        // Calculate the phi angle for this pixel
                        scalar pxRayPhi = iRayPhi - 0.5*deltaPhi + (j - 0.5)*deltaPhi/npPhi;

                        // Create a vector in the direction of this pixel
                        vector pixelDir = dom.anglesToDir(pxRayTheta, pxRayPhi);

                        // Get the cosine of the angle between the pixel direction and the surface normal
                        scalar cosA = pixelDir & surfNorm;

                        // Continue only if pixel ray is directed from other side of interface to this side
                        if (cosA > 0.0)
                        {
                            // Check that outgoing refracted ray is within the critical angle, otherwise
                            // no energy will be transferred across the interface
                            if ((1 - (1 - cosA*cosA)/(nRatio*nRatio)) > 0.0)
                            {
                                // Calculate direction and ID of refracted ray
                                vector refractDir = (pixelDir - cosA*surfNorm)/nRatio + Foam::sqrt(1 - (1 - cosA*cosA)/(nRatio*nRatio))*surfNorm;
                                label refractRay = dom.dirToRayId(refractDir, iBand);

                                // Accumulate the specular flux only if the outgoing ray matches this ray
                                if (refractRay == rayId)
                                {
                                    // Cosine of the refracted outgoing ray with
                                    // the surface normal -- exact continuous
                                    // value.
                                    scalar cosB = refractDir & surfNorm;

                                    // Calculate interface reflectivity
                                    scalar r1 = (nA*cosB - nB*cosA)/(nA*cosB + nB*cosA);
                                    scalar r2 = (nA*cosA - nB*cosB)/(nA*cosA + nB*cosB);
                                    scalar R = 0.5*(r1*r1 + r2*r2);

                                    // Calculate integral of direction vector wrt solid angle over pixel
                                    vector intDirOmega = dom.intDirOmega(pxRayTheta, pxRayPhi, deltaTheta/npTheta, deltaPhi/npPhi);

                                    // Accumulate the specular flux. Transmitted
                                    // radiance picks up an (nB/nA)^2 factor from
                                    // solid-angle compression across the interface.
                                    specular += NbrRaySet[iRay][faceI]*nRatio*nRatio*(1-R)*(intDirOmega & surfNorm);
                                }
                            }
                        }
                    }
                }
            }

            // Calculate the face value of the radiation intensity for this ray
            refValue()[faceI] = diffuseFraction_*diffuse/pi + (1.0 - diffuseFraction_)*specular/posIntDirOmega;
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = posIntDirOmega/magIntDirOmega;
        }
        // Case 2: Ray is going into wall with no overhang -> Treat as zero gradient
        else
        {
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0;
        }
    }

    UPstream::msgType() = oldTag;
    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::optical::multiBandTransInteriorCoupledFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("nBands") << nBands_ << token::END_STATEMENT << nl;
    os.writeKeyword("nNbg") << nNbg_ << token::END_STATEMENT << nl;
    os.writeKeyword("nOwn") << nOwn_ << token::END_STATEMENT << nl;
    os.writeKeyword("diffuseFraction") << diffuseFraction_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        makePatchTypeField
        (
            fvPatchScalarField,
            multiBandTransInteriorCoupledFvPatchScalarField
        );
    }
} // End namespace Foam

// ************************************************************************* //
