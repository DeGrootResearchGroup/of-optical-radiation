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

Class
    Foam::optical::radiationCoupledMixedFvPatchScalarField

Description
    Transparent radiation coupling for a multi-region interface where the
    two sides share the same refractive index. Replaces
    `refractiveCoupled` with `n_A = n_B` -- which still runs pixelation,
    candidate-ray gathering, Fresnel reflection and the n^2 etendue
    scaling even though all four collapse to no-ops at matched index.

    Behaviour at a face:
      * If the ray flows INTO this domain through the face (ray
        direction has positive component along the inward normal),
        the boundary value of `I` is fixed to the same ray's
        patch-internal value on the neighbour side, mapped via the
        usual `mappedPatchBase::fromNeighbour` machinery.
      * Otherwise (ray flows OUT), the BC is zero-gradient -- the
        DOM upwind discretisation uses the cell-internal value.

    No pixelation. No Fresnel. No n^2 scaling. Per-face complexity is
    O(1) instead of O(nPixelTheta * nPixelPhi * nAngle).

    Both sides must declare `n` (per-band) and the lists must agree
    within `nMatchTol` (default 1e-9 relative). A mismatch is a fatal
    error pointing the user at `refractiveCoupled` instead -- catches
    the common case of an interface that should refract but was wired
    up to the wrong BC.

    Example usage:
        myInterfacePatchName
        {
            type            radiationCoupled;
            nBands          2;
            n               (1.33 1.33);
            value           uniform 0;
        }

    Needs to be on an underlying mappedPatch (the same poly-patch
    family that `refractiveCoupled` uses; both BCs draw their
    cross-region data from `mappedPatchBase::fromNeighbour`).

SourceFiles
    radiationCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef optical_radiationCoupledMixedFvPatchScalarField_H
#define optical_radiationCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

namespace Foam
{
namespace optical
{

class radiationCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

        label nBands_;
        DynamicList<scalar, 1, 0> n_;

        //- Has the per-construction refractive-index sanity check
        //  run yet? Deferred to first updateCoeffs() because the
        //  neighbour patch field may not exist at our construction
        //  time (e.g. multi-region case-build order).
        mutable bool nValidated_;


    // Private helpers

        //- Look up the neighbour patch field of the same name and
        //  return true iff its `n` per-band list matches this BC's
        //  within nMatchTol relative. FatalError otherwise.
        void validateN_() const;


public:

    //- Runtime type information
    TypeName("radiationCoupled");

    //- Relative tolerance for the per-band refractive-index match.
    static constexpr scalar nMatchTol = 1.0e-9;


    // Constructors

        radiationCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        radiationCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        radiationCoupledMixedFvPatchScalarField
        (
            const radiationCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        radiationCoupledMixedFvPatchScalarField
        (
            const radiationCoupledMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new radiationCoupledMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Per-band refractive index on this side (read-only access).
        const DynamicList<scalar, 1, 0>& n() const { return n_; }

        virtual void updateCoeffs();

        virtual void write(Ostream&) const;
};

} // End namespace optical
} // End namespace Foam

#endif

// ************************************************************************* //
