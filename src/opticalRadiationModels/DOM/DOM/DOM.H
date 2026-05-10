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
    Foam::optical::DOM

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media.

    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)

    Angular convention (standard / ISO 80000-2):
        theta is measured from the +z axis, in [0, pi]
        phi   is measured from the +x axis (toward +y), in [0, 2*pi)
        Direction vector:
            s = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))

    Note: Murthy & Mathur (1998) measure phi from the +y axis instead,
    giving s = (sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta)).
    Their Eq. (14) integral formula for S_i has sin(phi) on the x-component
    and cos(phi) on the y-component as a consequence; the equivalent here
    in intDirOmega() has these swapped because the convention is rotated
    90 degrees about z. The integrals and ray directions are physically
    identical under either convention, but the iPhi index maps to a
    different physical direction in the two conventions.

SourceFiles
    DOM.C

\*---------------------------------------------------------------------------*/

#ifndef opticalRadiationModelDOM_H
#define opticalRadiationModelDOM_H

#include "radiationModel.H"
#include "ray.H"
#include "phaseFunctionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                           Class DOM Declaration
\*---------------------------------------------------------------------------*/

class DOM
:
    public radiationModel
{
    // Private data

        //- Irradiance field  [W/m^2]
        volScalarField G_;

        //- Number of solid angles in theta
        label nTheta_;

        //- Number of solid angles in phi
        label nPhi_ ;

        // - Number of solid angles
        label nAngle_;

        //- Total number of rays (1 per direction)
        label nRay_;

        //- Number of wavelength bands
        label nBand_;

        //- Number of pixels in phi
        label nPixelPhi_;

        //- Number of pixels in theta
        label nPixelTheta_;

        //- List of pointers to irradiance fields (one per band)
        PtrList<volScalarField> GLambda_;

        //- List of pointers to per-direction radiance rays
        PtrList<ray> IRay_;

        //- Per-ray radiance snapshot taken at the start of each outer
        //  iteration. The in-scatter source for ray i uses the snapshot
        //  values of all I_j (Jacobi update) rather than the
        //  partially-updated current I_j (Gauss-Seidel). The Jacobi
        //  form removes the order-dependent staleness in the scatter
        //  coupling that otherwise can drive an oscillation in the
        //  outer iteration on strongly-coupled cases (multi-band 3-D
        //  with anisotropic phase function).
        PtrList<volScalarField> ISnapshot_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        label maxIter_;

        //- phi interval
        scalar deltaPhi_;

        //- theta interface
        scalar deltaTheta_;

        //- Pointer to phase function model
        autoPtr<phaseFunctionModel> phaseFunctionModel_;

    // Private member functions

        DOM(const DOM&) = delete;
        void operator=(const DOM&) = delete;

        //- Check that dimension of the mesh is valid for simulation
        void checkDim_();

        //- Set up a radiance ray
        void setRay_
        (
            const label i,
            const label iBand,
            const label iAngle,
            const scalar theta,
            const scalar phi
        );

public:

    //- Runtime type information
    TypeName("DOM");

    // Constructors

        //- Construct from components
        DOM(const volScalarField& I);

    //- Destructor
    virtual ~DOM();


    // Member functions

        // Edit

            //- Solve opticalRadiation equation(s)
            void calculate();

            //- Read opticalRadiation properties dictionary
            bool read();

            //- Update G
            void updateG();

            //- Get the ray ID by decomposing a radiance field name
            label nameToRayId(const word& name) const;

            //- Get the theta angle from a direction vector
            scalar dirToTheta(const vector& dir) const;

            //- Get the phi angle from a direction vector
            scalar dirToPhi(const vector& dir) const;

            //- Get the direction vector from angles
            vector anglesToDir(const scalar& theta, const scalar& phi) const;

            //- Get the ray ID by decomposing a direction and band id
            label dirToRayId(const vector& dir, const label& iBand) const;

            //- Calculate integral of direction vector wrt solid angle
            //  (using specified deltaTheta and deltaPhi values)
            vector intDirOmega
            (
                const scalar& theta,
                const scalar& phi,
                const scalar& deltaTheta,
                const scalar& deltaPhi
            ) const;

            //- Calculate integral of direction vector wrt solid angle
            //  (using default deltaTheta and deltaPhi values)
            vector intDirOmega
            (
                const scalar& theta,
                const scalar& phi
            ) const;

            //- Ray radiance for rayI
            inline const ray& IRay(const label rayI) const;

            //- Number of angles in theta
            inline label nTheta() const;

            //- Number of angles in phi
            inline label nPhi() const;

            //- Number of angles in total
            inline label nAngle() const;

            //- Number of rays
            inline label nRay() const;

            //- Number of bands
            inline label nBand() const;

            //- Number of pixels in phi
            inline label nPixelPhi() const;

            //- Number of pixels in theta
            inline label nPixelTheta() const;

            //- Theta interval
            inline scalar deltaTheta() const;

            //- Phi interval
            inline scalar deltaPhi() const;

            //- Absorbtion coefficient
            inline const volScalarField& A(const label iLambda) const;

            //- Scattering coefficient
            inline const volScalarField& S(const label iLambda) const;

            //- Const access to irradiance field
            virtual const volScalarField& G() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
