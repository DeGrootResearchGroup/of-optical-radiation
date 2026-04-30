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
    Foam::optical::constantScatter

Description
    Constant opticalRadiation scatter coefficient

SourceFiles
    scatterModel.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicScattering_H
#define anisotropicScattering_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "scalar.H"

#include "inScatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                           Class constantScatter Declaration
\*---------------------------------------------------------------------------*/

class anisotropicScattering
:
    public inScatterModel
{
public:


    	DynamicList<scalar, 1, 0> phaseFuncCoef_;
	

private:

    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;
        
        //- Bands
        label nBands_;
        
        
        //- Bands
  //      FixedList<Vector2D<scalar>, maxBands_> iBands_;

/*
        //- Proportion of the heat released rate emitted
        FixedList<scalar, maxBands_> iEhrrCoeffs_;

*/

        // Total wave length covered by the bands
 //       scalar totalWaveLength_;



public:

    //- Runtime type information
    TypeName("anisotropicScattering");


    // Constructors

        //- Construct from components
        anisotropicScattering
        (
            const dictionary& dict
          //  , const fvMesh& mesh
        );


    // Destructor

       virtual ~anisotropicScattering();


    // Member Operators


/* 
        //- Lower and upper limit of band i
        inline const Vector2D<scalar>& bands(const label i) const
        {
            return iBands_[i];
        }
            
       inline bool isGrey() const
        {
            return false;
        }
        * 
        //- Add contribution of ILambda to the total radiative intensity in
        //  direction i
        tmp<volScalarField> addIntensity
        (
            const label i,
            const volScalarField& ILambda
        ) const;
*/
         Foam::scalar correct
        (
          const  scalar  rayCos,
            const label  iBand
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
