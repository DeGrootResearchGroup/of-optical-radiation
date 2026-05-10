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
    Foam::optical::mieModel

Description
    Full Mie phase function for monodisperse spherical particles. Per band,
    constructs a Mie kernel from the size parameter and complex relative
    refractive index. The phaseShape returns the raw BHMIE intensity
    |S_1|^2 + |S_2|^2 (unnormalised); the row-normalisation in
    buildPhaseTable absorbs the absolute scaling.

    Reads the same physical parameters as mieExtinction (radius, mParticle,
    mMedium, wavelengths). Keeping them duplicated rather than coupled
    matches how every other (extinctionModel, phaseFunctionModel) pair is
    configured -- the two objects remain independently selectable.

SourceFiles
    mieModel.C

\*---------------------------------------------------------------------------*/

#ifndef mieModel_H
#define mieModel_H

#include "phaseFunctionModel.H"
#include "mieKernel.H"
#include "PtrList.H"

#include <complex>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

class DOM;

/*---------------------------------------------------------------------------*\
                          Class mieModel Declaration
\*---------------------------------------------------------------------------*/

class mieModel
:
    public phaseFunctionModel
{
private:

    // Private data

        //- Particle radius [m]
        scalar radius_;

        //- Complex particle refractive index
        std::complex<scalar> mParticle_;

        //- Real refractive index of the host medium
        scalar mMedium_;

        //- Per-band central wavelengths in vacuum [nm]
        List<scalar> wavelengths_;

        //- One Mie kernel per band; constructed once in the ctor and
        //  consulted by phaseShape() during table build.
        PtrList<mieKernel> kernels_;


    // Private member functions

        //- Override of phaseFunctionModel::phaseShape: returns the raw
        //  BHMIE intensity |S_1|^2 + |S_2|^2 evaluated at cos theta.
        //  Absolute scaling does not matter (row-normalisation absorbs
        //  it); only the angular shape contributes to the in-scatter
        //  source.
        virtual scalar phaseShape
        (
            const scalar cosV,
            const label iBand
        ) const override;


public:

    //- Runtime type information
    TypeName("mieModel");


    // Constructors

        mieModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~mieModel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
