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
    Foam::solarSurfaceMixedFvPatchScalarField

Description
    opticalRadiation temperature specified

SourceFiles
    solarSurfaceMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef solarSurfaceMixedFvPatchScalarField_H
#define solarSurfaceMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{
/*---------------------------------------------------------------------------*\
        Class solarSurfaceMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solarSurfaceMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- incident opticalRadiation angle distribution
        scalar  nNbg_ ;
        scalar  nOwn_; 
        scalar  diffuseCoeff_;
         
        scalar  directBeam_;              // incident opticalRadiation intensity
        scalar  skyDiffuse_;   
        scalar  groundReflected_;  
        scalar  zenithAngle_;         // start from the vertical 
        scalar  azimuthAngle_;      // start from the east
        
		label   nBands_;  
    	DynamicList<scalar,1,0>  opticalRadiationBandDist_;
    	
//		autoPtr<scalarField>   refraction;
	//	autoPtr<scalarField>   diffusiveRefraction;		//	        DynamicList<scalar,1,0>  angleDist; 
//		label  initFlag_ ; 		

		        
public:

    //- Runtime type information
    TypeName("solarSurface");


    // Constructors

        //- Construct from patch and internal field
        solarSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solarSurfaceMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  solarSurfaceMixedFvPatchScalarField onto a new patch
        solarSurfaceMixedFvPatchScalarField
        (
            const solarSurfaceMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        solarSurfaceMixedFvPatchScalarField
        (
            const solarSurfaceMixedFvPatchScalarField&,
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
                new solarSurfaceMixedFvPatchScalarField(*this, iF)
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
