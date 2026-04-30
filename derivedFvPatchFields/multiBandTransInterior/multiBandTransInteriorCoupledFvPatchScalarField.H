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

Class
    multiBandTransInteriorCoupledFvPatchScalarField

Description
    Mixed boundary condition for radiation at a transmissive interior interface.

    Example usage:
        myInterfacePatchName
        {
            type            multiBandTransInteriorCoupled;
            diffuseFraction 0.0;
            nBands          2;
            nNbg            (1.5 1.5);
            nOwn            (1.0 1.0);
            value           $internalField;
        }

    Needs to be on underlying mappedPatch.

SourceFiles
    multiBandTransInteriorCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef multiBandTransInteriorCoupledFvPatchScalarField_H
#define multiBandTransInteriorCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {

        class multiBandTransInteriorCoupledFvPatchScalarField
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
            TypeName("multiBandTransInteriorCoupled");


            // Constructors

            //- Construct from patch and internal field
            multiBandTransInteriorCoupledFvPatchScalarField
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&
            );

            //- Construct from patch, internal field and dictionary
            multiBandTransInteriorCoupledFvPatchScalarField
            (
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const dictionary&
            );

            //- Construct by mapping given
            //  multiBandTransInteriorCoupledFvPatchScalarField onto a new patch
            multiBandTransInteriorCoupledFvPatchScalarField
            (
                const multiBandTransInteriorCoupledFvPatchScalarField&,
                const fvPatch&,
                const DimensionedField<scalar, volMesh>&,
                const fvPatchFieldMapper&
            );

            //- Construct as copy setting internal field reference
            multiBandTransInteriorCoupledFvPatchScalarField
            (
                const multiBandTransInteriorCoupledFvPatchScalarField&,
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
                    new multiBandTransInteriorCoupledFvPatchScalarField
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

    } // End namespace photoBio

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
