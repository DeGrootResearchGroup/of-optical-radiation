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
    Foam::photoBio

Description
    Namespace for photoBio modelling

Class
    Foam::photoBio::photoBioModel

Description
    Top level model for photoBio modelling

SourceFiles
    photoBioModel.C

\*---------------------------------------------------------------------------*/

#ifndef photoBioModel_H
#define photoBioModel_H

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
namespace photoBio
{

// Forward declaration of classes
class extinctionModel;
class inScatterModel;

/*---------------------------------------------------------------------------*\
                       Class photoBioModel Declaration
\*---------------------------------------------------------------------------*/

class photoBioModel
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
        Switch photoBio_;

        //- photoBio model dictionary
        dictionary coeffs_;

        //- References to the photoBio sub-models
         autoPtr<extinctionModel> extinction_;
         autoPtr<inScatterModel> inScatter_;

private:

    // Private Member Functions

        photoBioModel(const photoBioModel&) = delete;
        void operator=(const photoBioModel&) = delete;

public:

    //- Runtime type information
    TypeName("photoBioModel");

    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             photoBioModel,
             dictionary,
             (
                 const volScalarField& intensity
             ),
             (intensity)
         );


    // Constructors

        //- Null constructor
        photoBioModel(const volScalarField& intensity);

        //- Construct from components
        photoBioModel(const word& type, const volScalarField& intensity);


    // Selectors

         //- Return a reference to the selected photoBio model
        static autoPtr<photoBioModel> New(const volScalarField& intensity);


    //- Destructor
    virtual ~photoBioModel();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

            //- Solve photoBio equation(s)
            virtual void calculate() = 0;

            //- Read photoBioProperties dictionary
            virtual bool read() = 0;

            virtual const volScalarField& G() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
