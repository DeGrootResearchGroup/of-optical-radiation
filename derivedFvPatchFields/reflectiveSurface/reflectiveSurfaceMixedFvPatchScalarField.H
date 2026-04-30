/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::semiTransSurfaceMixedFvPatchScalarField

Description
    photoBio temperature specified

SourceFiles
    semiTransSurfaceMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef reflectiveSurfaceMixedFvPatchScalarField_H
#define reflectiveSurfaceMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{
/*---------------------------------------------------------------------------*\
        Class reflectiveSurfaceMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class reflectiveSurfaceMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

		scalar  diffuseFraction_; 
        scalar  reflectionCoef_; 
   //     DynamicList<scalar,1,0>  angleDist; 
        
        
public:

    //- Runtime type information
    TypeName("reflectiveSurface");


    // Constructors

        //- Construct from patch and internal field
        reflectiveSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        reflectiveSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  reflectiveSurfaceMixedFvPatchScalarField onto a new patch
        reflectiveSurfaceMixedFvPatchScalarField
        (
            const reflectiveSurfaceMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        reflectiveSurfaceMixedFvPatchScalarField
        (
            const reflectiveSurfaceMixedFvPatchScalarField&,
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
                new reflectiveSurfaceMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access
             
        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // I-O

            //- Write
             void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
