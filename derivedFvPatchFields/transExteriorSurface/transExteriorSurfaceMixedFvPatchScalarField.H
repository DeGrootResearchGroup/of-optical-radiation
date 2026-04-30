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
    Foam::transExteriorSurfaceMixedFvPatchScalarField

Description
    photoBio temperature specified

SourceFiles
    transExteriorSurfaceMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef transExteriorSurfaceMixedFvPatchScalarField_H
#define transExteriorSurfaceMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{
/*---------------------------------------------------------------------------*\
        Class transExteriorSurfaceMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class transExteriorSurfaceMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- incident photoBio angle distribution
        scalar  nNbg_ ;
        scalar  nOwn_;       
        scalar  I0_;              // incident photoBio intensity
        label   nBands_;
        scalar  diffuseFraction_;    
        scalar  beamWidthPhi_	;
        scalar  beamWidthTheta_	;
        vector  beamDir_ ;         // incident photoBio direction
        bool beamNormToSurf_;
		label  initFlag_ ; 
    	DynamicList<scalar,1,0>  photoBioBandDist_;
		autoPtr<scalarField>   specularRefraction;
		autoPtr<scalarField>   diffusiveRefraction;	
		
	//	        DynamicList<scalar,1,0>  angleDist; 
		        
		scalar  beamAnglePhi_;
		scalar  beamAngleTheta_;
		
public:

    //- Runtime type information
    TypeName("transExteriorSurface");


    // Constructors

        //- Construct from patch and internal field
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  transExteriorSurfaceMixedFvPatchScalarField onto a new patch
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const transExteriorSurfaceMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        transExteriorSurfaceMixedFvPatchScalarField
        (
            const transExteriorSurfaceMixedFvPatchScalarField&,
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
                new transExteriorSurfaceMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
			virtual void dirToAngle(const vector& ,scalar&, scalar& ) const;
			
        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
