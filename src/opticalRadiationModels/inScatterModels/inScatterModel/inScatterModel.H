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
    Foam::optical::scatterModel

Description
    Base class for opticalRadiation scattering

\*---------------------------------------------------------------------------*/

#ifndef inScatterModel_H
#define inScatterModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                           Class inScatterModel Declaration
\*---------------------------------------------------------------------------*/

class inScatterModel
{

protected:

    // Protected data

        //- Reference to the fvMesh
   //     const fvMesh& mesh_;
        
        
        //- opticalRadiation model dictionary
        const dictionary dict_;
        
        
public:

    //- Runtime type information
    TypeName("inScatterModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            inScatterModel,
            dictionary,
            (
                const dictionary& dict
  // ,            const fvMesh& mesh
            ),
            (dict)
   //         (dict, mesh)
        );


    // Constructors

        //- Construct from components
        inScatterModel
        (
            const dictionary& dict
    //        , const fvMesh& mesh
        );


    // Selector

        static autoPtr<inScatterModel> New
        (
            const dictionary& dict
    //        ,  const fvMesh& mesh
        );


    // Destructor

        virtual ~inScatterModel();


    // Member Functions

        // Access
            //- Reference to the mesh
  /*          inline const fvMesh& mesh() const
            {
                return mesh_;
            }
*/
            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }
        
            virtual bool inScatter() const;
                        
/*
            virtual scalar correct
            (
				scalar  angle,
				label  iBand
            ) const;
            */
            
            virtual void update
            (
				label  nPhi,
				label  nTheta,
				label  nBand,
				label  nDim
            ) const;
            
            
            virtual scalar correct
            (
				label  rayI,
				label  rayJ,
				label  iBand
            ) const;


                
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

