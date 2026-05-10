/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::optical::refractiveCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for radiation at a refractive interior
    interface (mapped patch). Splits arriving radiance into Fresnel
    transmission to the neighbour region and reflection back into the
    own region, including the étendue n² factor.

    Example usage:
        myInterfacePatchName
        {
            type            refractiveCoupled;
            diffuseFraction 0.0;
            nBands          2;
            nNbg            (1.5 1.5);
            nOwn            (1.0 1.0);
            value           $internalField;
        }

    Needs to be on underlying mappedPatch.

SourceFiles
    refractiveCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef optical_refractiveCoupledMixedFvPatchScalarField_H
#define optical_refractiveCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {

        class refractiveCoupledMixedFvPatchScalarField
        :
        public mixedFvPatchScalarField
        {

            // Private data

            scalar diffuseFraction_;
            label nBands_;
            DynamicList<scalar,1,0> nNbg_;
            DynamicList<scalar,1,0> nOwn_;

        public:

            //- Runtime type information
            TypeName("refractiveCoupled");


            // Constructors

            //- Construct from patch and internal field
            refractiveCoupledMixedFvPatchScalarField
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&
            );

            //- Construct from patch, internal field and dictionary
            refractiveCoupledMixedFvPatchScalarField
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const dictionary&
            );

            //- Construct by mapping given
            //  refractiveCoupledMixedFvPatchScalarField onto a new patch
            refractiveCoupledMixedFvPatchScalarField
            (
                const refractiveCoupledMixedFvPatchScalarField&,
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const fvPatchFieldMapper&
            );

            //- Construct as copy setting internal field reference
            refractiveCoupledMixedFvPatchScalarField
            (
                const refractiveCoupledMixedFvPatchScalarField&,
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
                    new refractiveCoupledMixedFvPatchScalarField
                    (
                        *this,
                        iF
                    )
                );
            }


            // Member functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
        };

    } // End namespace optical

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
