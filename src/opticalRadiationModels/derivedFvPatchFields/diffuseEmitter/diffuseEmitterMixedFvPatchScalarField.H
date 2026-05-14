/*---------------------------------------------------------------------------*\
    =========                 |
    \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
     \\    /   O peration     |
      \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
       \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
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

Class
    Foam::optical::diffuseEmitterMixedFvPatchScalarField

Description
    Boundary condition for a Lambertian (diffuse) emitter. Each face on
    the patch emits radiation with a direction-independent radiance into
    every ray bin that points away from the wall into the domain. Faces
    whose surface normal sees the ray direction pointing back into the
    wall get zero-gradient treatment, so the BC does not interact with
    radiation arriving at the patch.

    Mathematical model. For each ray d in the angular grid:
    \verbatim
        cos_alpha = (-n_face) . d           // n_face points out of domain
        if cos_alpha > 0:                   // ray points into domain
            L_d = emissivePower[iBand] / pi
        else:                               // ray points out of domain
            zero-gradient
    \endverbatim
    The 1/pi factor is the Lambertian normalisation
    (integral over the outward hemisphere of cos(theta) dOmega = pi), so
    `emissivePower` is the radiometric flux [W/m^2] emitted by the
    patch and `L_d` is the resulting radiance [W/m^2/sr] -- units match
    the radiative-transfer equation that DOM is solving.

    Dictionary entries:
    \verbatim
        radSource
        {
            type            diffuseEmitter;
            nBands          2;
            emissivePower   (5 5);          // W/m^2 per band
            value           uniform 0;
        }
    \endverbatim

    Conventions:
      - `emissivePower` has length `nBands`. Per-band entries let a
        single emitter feed a multi-band DOM solve with different
        intensities in each band (e.g. a UV-C lamp tabulated by
        emission line).
      - The patch is treated as a planar Lambertian aperture: the
        emitted radiance is constant across direction and uniform over
        every face. For a non-uniform emissivePower (e.g. spatially
        varying lamp output across a window), pre-split the patch into
        multiple sub-patches with separate dictionary entries.
      - This BC has no reflective component: radiation arriving at the
        patch from the domain is absorbed (zero-gradient on the
        into-wall hemisphere, set by `valueFraction = 0`). For a
        purely reflecting wall (no emission), use the `reflective` BC.
        Combined emit+reflect surfaces are not currently supported in
        a single BC instance.

SourceFiles
    diffuseEmitterMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef optical_diffuseEmitterMixedFvPatchScalarField_H
#define optical_diffuseEmitterMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
  Class diffuseEmitterMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class diffuseEmitterMixedFvPatchScalarField : public mixedFvPatchScalarField
{
    // Private data

    label nBands_;
    DynamicList<scalar,1,0>  emissivePower_;

public:

    //- Runtime type information
    TypeName("diffuseEmitter");

    // Constructors

    //- Construct from patch and internal field
    diffuseEmitterMixedFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    diffuseEmitterMixedFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given a diffuseEmitterMixedFvPatchScalarField onto a new patch
    diffuseEmitterMixedFvPatchScalarField
    (
        const diffuseEmitterMixedFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy setting internal field reference
    diffuseEmitterMixedFvPatchScalarField
    (
        const diffuseEmitterMixedFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return tmp<fvPatchScalarField>
        (
            new diffuseEmitterMixedFvPatchScalarField(*this, iF)
        );
    }

    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
