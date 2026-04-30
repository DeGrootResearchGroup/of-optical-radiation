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
    transInteriorSurfaceCoupledFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used by the
    conjugate heat transfer solver.
    Both sides use a mix of zeroGradient and neighbourvalue.

    Example usage:
        myInterfacePatchName
        {
            type                transInteriorSurfaceCoupled;
            neighbourFieldName  T;

        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

SourceFiles
    transInteriorSurfaceCoupledFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef transInteriorSurfaceCoupledFvPatchScalarField_H
#define transInteriorSurfaceCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{
/*---------------------------------------------------------------------------*\
        Class transInteriorSurfaceCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class transInteriorSurfaceCoupledFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region

        
        scalar  nNbg_;          // index of refraction 
        scalar  nOwn_;   
        scalar  diffuseFraction_;     
  //      DynamicList<scalar,1,0>  angleDist; 
        
public:

    //- Runtime type information
    TypeName("transInteriorSurfaceCoupled");


    // Constructors

        //- Construct from patch and internal field
        transInteriorSurfaceCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transInteriorSurfaceCoupledFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  transInteriorSurfaceCoupledFvPatchScalarField onto a new patch
        transInteriorSurfaceCoupledFvPatchScalarField
        (
            const transInteriorSurfaceCoupledFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        transInteriorSurfaceCoupledFvPatchScalarField
        (
            const transInteriorSurfaceCoupledFvPatchScalarField&,
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
                new transInteriorSurfaceCoupledFvPatchScalarField
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace photoBio

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
