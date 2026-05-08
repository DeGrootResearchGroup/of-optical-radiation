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
    Foam::optical::ray

Description
    Radiance for a ray in a given direction (units W/m^2/sr).

SourceFiles
    ray.C

\*---------------------------------------------------------------------------*/

#ifndef ray_H
#define ray_H

#include "surfaceFields.H"

#include "extinctionModel.H"
#include "phaseFunctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

// Forward declaration of classes
class DOM;

/*---------------------------------------------------------------------------*\
                    Class ray Declaration
\*---------------------------------------------------------------------------*/

class ray
{
public:

    static const word namePrefix;

private:

    // Private data

        //- Refence to the owner DOM object
        const DOM& dom_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Ray radiance [W/m^2/sr]
        autoPtr<volScalarField>  I_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Wavelength band index
        label iBand_;

        //- Angle index
        label iAngle_;

    // Private member functions

        ray(const ray&) = delete;
        void operator=(const ray&) = delete;

        //- Fill the upwind (Ji0) and downwind (Ji1) face-flux fields
        //  by summing contributions from each pixel of the ray's
        //  solid angle. Ji0 collects pixels with positive face-normal
        //  component (face is "in" for the pixel direction) and Ji1
        //  the negative ones; faces exactly tangent to the pixel
        //  direction (dot product zero) are dropped from both -- the
        //  Murthy-Mathur tangent-pixel choice, second-order accurate
        //  in the pixel discretisation. Both fields are filled in a
        //  single pixel-loop pass with no intermediate
        //  surfaceScalarField allocations.
        void computeFluxCoeffs_
        (
            surfaceScalarField& Ji0,
            surfaceScalarField& Ji1
        ) const;

public:

    // Constructors

        //- Construct from components
        ray
        (
            const DOM& dom,
            const fvMesh& mesh,
            const label  iBand,
            const label  iAngle,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta

        );

    // Destructor
    ~ray();

    // Member functions

        // Edit

            //- Update radiance on i direction
            scalar correct();

            //- Update boundaries
            void updateBoundary();

        // Access

            //- Return radiance
             inline const volScalarField& I() const;

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

             //- Return the index of band
            inline  label iBand() const;

            //- Return the index of the angle
            inline label iAngle() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
