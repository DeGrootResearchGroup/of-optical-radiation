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
    Foam::optical::iesEmitterMixedFvPatchScalarField

Description
    Boundary condition representing a real luminaire whose far-field
    angular intensity distribution is described by an IESNA LM-63
    Type C photometric file. Targets KrCl excimer fixtures, but
    works for any directional emitter packaged with an IES file.

    The IES table supplies the angular SHAPE of the emission. Absolute
    magnitude is provided separately as a per-band total radiant flux
    P [W] -- this lets the same BC consume IES files in either
    photometric (candela) or radiometric (W/sr) units, since the
    table is renormalised against P at runtime.

    Mathematical model. For each ray d in the angular grid going INTO
    the domain through the patch:
    \verbatim
        L_d = (P / (A_patch * Phi_table)) * I_table(d) / max(d.n_avg, eps)
    \endverbatim
    where:
      - A_patch  = global patch area (summed across processors)
      - n_avg    = patch-averaged inward normal
                   = -sum(Sf) / sum(magSf) (inward into the domain)
      - Phi_table = sum over outgoing rays of I_table(d) * Omega_d
                    (no cos weighting -- normalisation by intensity, not
                     radiance)
      - eps       = small floor on cos(d, n_avg) to bound L_d at grazing
    This choice gives total emitted radiometric flux exactly equal to
    P, and per-direction emitted intensity I_emit(d) = L_d * A_patch *
    cos(d, n_avg) proportional to I_table(d) -- i.e. the patch's
    far-field intensity pattern matches the IES file shape. See
    CLAUDE.md "iesEmitter" for the derivation.

    Outgoing rays at the patch (cos(d, n_face) <= 0) get zero-gradient
    treatment, identical to the other emitter BCs.

    Dictionary entries:
    \verbatim
        louver
        {
            type            iesEmitter;
            nBands          1;
            iesFile         "constant/krcl-louvered.ies";
            power           (12.5);                   // W per band
            fixtureAxis     (0 0 -1);                 // IES gamma=0 dir
            fixtureUp       (1 0 0);                  // IES h=0 dir
            value           uniform 0;
        }
    \endverbatim

    Conventions:
      - fixtureAxis is the direction in global (x,y,z) coordinates
        corresponding to IES vertical angle gamma = 0 (the fixture's
        nominal beam axis). For a downlight pointing -z, this is
        (0, 0, -1). The BC normalises it.
      - fixtureUp is the direction in global coordinates corresponding
        to IES horizontal angle h = 0. It is orthogonalised against
        fixtureAxis at construction. For axisymmetric IES files
        (single horizontal angle) the choice does not matter; supply
        any vector not parallel to fixtureAxis.
      - The patch is treated as a planar emitter aperture. A non-planar
        patch is handled (using the patch-averaged normal for the
        radiance normalisation), but the BC's intent is for nominally
        flat luminaire windows.
      - Multi-band: the same IES shape is applied across all bands,
        scaled by the per-band entry of `power`. If a fixture has
        wavelength-dependent angular structure, declare the patch
        once per band with separate IES files (not yet supported in
        a single BC instance).

SourceFiles
    iesEmitterMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef iesEmitterMixedFvPatchScalarField_H
#define iesEmitterMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "iesPhotometry.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
       Class iesEmitterMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class iesEmitterMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Path to the IES Type C photometric file
        fileName iesFile_;

        //- Number of wavelength bands
        label nBands_;

        //- Total radiant flux per band emitted into the domain [W]
        DynamicList<scalar, 1, 0> power_;

        //- Direction (global frame) of IES vertical angle gamma = 0.
        //  Normalised at construction.
        vector fixtureAxis_;

        //- Direction (global frame) of IES horizontal angle h = 0.
        //  Orthogonalised against fixtureAxis_ and normalised at
        //  construction.
        vector fixtureUp_;


    // Lazy-initialised cache (filled on first updateCoeffs())

        //- Parsed IES table
        mutable autoPtr<iesPhotometry> ies_;

        //- True once cachedL_ has been computed for this BC's ray
        mutable bool cacheValid_;

        //- Radiance assigned to this BC's ray on outgoing-into-domain
        //  faces [W/m^2/sr]. Zero if this ray's direction is not
        //  emitted by the fixture.
        mutable scalar cachedL_;


    // Private member functions

        //- Compute the IES vertical angle gamma [degrees] for a global
        //  direction d
        scalar gammaDegFromDir_(const vector& d) const;

        //- Compute the IES horizontal angle h [degrees] for a global
        //  direction d
        scalar hDegFromDir_(const vector& d) const;


public:

    //- Runtime type information
    TypeName("iesEmitter");


    // Constructors

        //- Construct from patch and internal field
        iesEmitterMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field, and dictionary
        iesEmitterMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        iesEmitterMixedFvPatchScalarField
        (
            const iesEmitterMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        iesEmitterMixedFvPatchScalarField
        (
            const iesEmitterMixedFvPatchScalarField&,
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
                new iesEmitterMixedFvPatchScalarField(*this, iF)
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
