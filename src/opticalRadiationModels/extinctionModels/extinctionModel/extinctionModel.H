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
    Foam::optical::extinctionModel

Description
    Model to supply absorption and scattering coefficients

\*---------------------------------------------------------------------------*/

#ifndef optical_extinctionModel_H
#define optical_extinctionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                   Class extinctionModel Declaration
\*---------------------------------------------------------------------------*/

class extinctionModel
{

protected:

    // Protected data

        //- Model dictionary
        const dictionary dict_;

        //- Reference to the fvMesh
        const fvMesh& mesh_;

        //- Number of bands
        label nBands_;

        //- List of pointers to absorption coefficient fields (one per band)
        PtrList<volScalarField> ALambda_;

        //- List of pointers to scattering coefficient fields (one per band)
        PtrList<volScalarField> SLambda_;

public:

    //- Runtime type information
    TypeName("extinctionModel");

    //- Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            extinctionModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );

    // Constructors

        //- Construct from components
        extinctionModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Selector
    static autoPtr<extinctionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~extinctionModel();

    // Member Functions

        // Initialization

            //- Initialize extinction coefficients based on number of bands
            void init(const label nBands);

        // Access

            //- Reference to the mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Reference to the dictionary
            inline const dictionary& dict() const
            {
                return dict_;
            }

            //- Number of bands
            inline label nBands() const
            {
                return nBands_;
            }

            //- Absorbtion coefficient
            inline const volScalarField& A(const label iLambda) const
            {
                return ALambda_[iLambda];
            }

            //- Scattering coefficient
            inline const volScalarField& S(const label iLambda) const
            {
                return SLambda_[iLambda];
            }

        // Edit

            //- Correct extinction coefficients
            virtual void correct() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

