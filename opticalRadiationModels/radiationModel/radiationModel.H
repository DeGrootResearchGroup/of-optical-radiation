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

Namespace
    Foam::optical

Description
    Namespace for opticalRadiation modelling

Class
    Foam::optical::radiationModel

Description
    Top level model for opticalRadiation modelling

SourceFiles
    radiationModel.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModel_H
#define radiationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "extinctionModel.H"
#include "inScatterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

// Forward declaration of classes
class extinctionModel;
class inScatterModel;

/*---------------------------------------------------------------------------*\
                       Class radiationModel Declaration
\*---------------------------------------------------------------------------*/

class radiationModel
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Model specific dictionary input parameters
        Switch opticalRadiation_;

        //- opticalRadiation model dictionary
        dictionary coeffs_;

        //- References to the opticalRadiation sub-models
         autoPtr<extinctionModel> extinction_;
         autoPtr<inScatterModel> inScatter_;

private:

    // Private Member Functions

        radiationModel(const radiationModel&) = delete;
        void operator=(const radiationModel&) = delete;

public:

    //- Runtime type information
    TypeName("radiationModel");

    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             radiationModel,
             dictionary,
             (
                 const volScalarField& intensity
             ),
             (intensity)
         );


    // Constructors

        //- Null constructor
        radiationModel(const volScalarField& intensity);

        //- Construct from components
        radiationModel(const word& type, const volScalarField& intensity);


    // Selectors

         //- Return a reference to the selected opticalRadiation model
        static autoPtr<radiationModel> New(const volScalarField& intensity);


    //- Destructor
    virtual ~radiationModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve opticalRadiation equation(s)
            virtual void calculate() = 0;

            //- Read opticalRadiationProperties dictionary
            virtual bool read() = 0;

            virtual const volScalarField& G() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
