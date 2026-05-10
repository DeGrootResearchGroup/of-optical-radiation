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

Class
    Foam::optical::collimatedBeamMixedFvPatchScalarField

Description
    Boundary condition for a collimated radiation beam entering the domain
    through a wall. Approximated within the discrete-ordinates angular
    discretisation by assigning the beam radiance to the single ray bin
    that contains the beam direction; all other rays at the patch see
    zero incoming radiance (they are absorbing for radiation arriving
    along their direction).

    Inspired by the irradiation-beam treatment in ANSYS Fluent's DO
    model (Theory Guide §5.3.6.7.5). The current form is the simplest
    delta-direction beam: one discrete ray carries the beam, no spread
    over neighbouring rays.

    The beam strength may be specified either directly as a radiance
    or via a per-area irradiation flux:

    Mode 1 -- direct radiance:
    \verbatim
        sourceWall
        {
            type            collimatedBeam;
            nBands          1;
            beamDirection   (0.9239 0.3827 0);
            beamRadiance    (10);            // W/m^2/sr per band
            value           uniform 0;
        }
    \endverbatim

    Mode 2 -- irradiation flux (Fluent-style):
    \verbatim
        sourceWall
        {
            type            collimatedBeam;
            nBands          1;
            beamDirection   (0.9239 0.3827 0);
            emissivePower   (1000);          // W/m^2 per band, beam-normal
            value           uniform 0;
        }
    \endverbatim

    Exactly one of `beamRadiance` and `emissivePower` must be set.

    Conventions:
    - beamDirection is interpreted in the global (x,y,z) frame and is
      the direction the beam *travels* (i.e. into the domain through
      the wall). The BC normalises it.
    - In mode 1, beamRadiance is the radiance L (W/m^2/sr) assigned to
      the discrete ray containing beamDirection.
    - In mode 2, emissivePower is the per-area flux carried by the
      beam, measured through a surface *perpendicular to
      beamDirection* (matches Fluent's "Apply Direct Irradiation
      Parallel to the Beam" convention). The BC converts to a per-ray
      radiance via L = emissivePower / omega_ray, where omega_ray is
      the chosen discrete ray's solid angle. The result is independent
      of face orientation, which is consistent with a true collimated
      beam.
    - The BC implicitly assumes the beam is narrow enough to fit in a
      single discrete-ordinates bin. Within the DOM's angular
      resolution this is the only regime that's representable;
      spreading energy across multiple rays for a wider beam is not
      yet implemented.
    - For faces where beamDirection points away from the domain
      (cos(beamDir, surfNorm_into_domain) <= 0), the beam ray is set
      to zero -- i.e. the beam only illuminates faces it actually
      points into.

SourceFiles
    collimatedBeamMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef optical_collimatedBeamMixedFvPatchScalarField_H
#define optical_collimatedBeamMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
       Class collimatedBeamMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class collimatedBeamMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Direction the beam travels (unit vector; normalised at read)
        vector beamDirection_;

        //- Number of wavelength bands
        label nBands_;

        //- Radiance per band (W/m^2/sr) carried by the discrete ray
        //  bin containing beamDirection_. Set in mode 1; size 0 in mode 2.
        DynamicList<scalar, 1, 0> beamRadiance_;

        //- Beam-normal flux per band (W/m^2). Set in mode 2; size 0 in
        //  mode 1. Converted to a per-ray radiance via
        //  L = emissivePower / omega_ray at updateCoeffs time.
        DynamicList<scalar, 1, 0> emissivePower_;


public:

    //- Runtime type information
    TypeName("collimatedBeam");


    // Constructors

        //- Construct from patch and internal field
        collimatedBeamMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field, and dictionary
        collimatedBeamMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        collimatedBeamMixedFvPatchScalarField
        (
            const collimatedBeamMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        collimatedBeamMixedFvPatchScalarField
        (
            const collimatedBeamMixedFvPatchScalarField&,
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
                new collimatedBeamMixedFvPatchScalarField(*this, iF)
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
