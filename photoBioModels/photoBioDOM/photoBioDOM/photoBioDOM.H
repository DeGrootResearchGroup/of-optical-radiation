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
    Foam::photoBio::photoBioDOM

Description

    Finite Volume Discrete Ordinates Method. Solves the RTE equation for n
    directions in a participating media.

    In 1D the direction of the rays is X (nPhi and nTheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only nPhi is considered)
    In 3D (nPhi and nTheta are considered)

SourceFiles
    photoBioDOM.C

\*---------------------------------------------------------------------------*/

#ifndef photoBioModelphotoBioDOM_H
#define photoBioModelphotoBioDOM_H

#include "photoBioModel.H"
#include "photoBioIntensityRay.H"
#include "phaseFunctionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace photoBio
{

/*---------------------------------------------------------------------------*\
                           Class photoBioDOM Declaration
\*---------------------------------------------------------------------------*/

class photoBioDOM
:
    public photoBioModel
{
    // Private data

        //- Irradiance field  [W/m^2]
        volScalarField G_;

        //- Inscattering field
        volScalarField diffusionScatter_;

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

        //- List of phase function values (one per band)
        DynamicList<scalar, 1, 0> pf0_;

        //- List of pointers to intensity rays
        PtrList<photoBioIntensityRay> IRay_;

        //- Convergence criterion
        scalar convergence_;

        //- Maximum number of iterations
        scalar maxIter_;

        //- phi interval
        scalar deltaPhi_;

        //- theta interface
        scalar deltaTheta_;

        //- Pointer to phase function model
        autoPtr<phaseFunctionModel> phaseFunctionModel_;

    // Private member functions

        photoBioDOM(const photoBioDOM&) = delete;
        void operator=(const photoBioDOM&) = delete;

        //- Check that dimension of the mesh is valid for simulation
        void checkDim_();

        //- Set up an intensity ray
        void setRay_(const label i,
                     const label iBand,
                     const label iAngle,
                     const scalar phi,
                     const scalar theta);

public:

    //- Runtime type information
    TypeName("photoBioDOM");

    // Constructors

        //- Construct from components
        photoBioDOM(const volScalarField& intensity);

    //- Destructor
    virtual ~photoBioDOM();


    // Member functions

        // Edit

            //- Solve photoBio equation(s)
            void calculate();

            //- Read photoBio properties dictionary
            bool read();

            //- Update G
            void updateG();

            //- Get the ray ID by decomposing an intensity field name
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

            //- Ray intensity for rayI
            inline const photoBioIntensityRay& IRay(const label rayI) const;

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

            //- Phase function values
            inline scalar pf0(const label iLambda) const;

            //- Const access to irradiance field
            virtual const volScalarField& G() const;

            //- Const access to
            inline const volScalarField& diffusionScatter() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "photoBioDOMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace photoBio
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
