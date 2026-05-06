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
    refractive index, then precomputes the row-normalised pixel-integrated
    P[i,j,band] table on the same sub-angle scheme as
    HenyeyGreensteinModel / rayleighModel.

    Reads the same physical parameters as mieExtinction (radius, mParticle,
    mMedium, wavelengths). Keeping them duplicated rather than coupled
    matches how every other (extinctionModel, phaseFunctionModel) pair is
    configured -- the two objects remain independently selectable.

SourceFiles
    mieModel.C

\*---------------------------------------------------------------------------*/

#ifndef mieModel_H
#define mieModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"

#include "phaseFunctionModel.H"

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

class mieModel : public phaseFunctionModel
{
private:

    // Private data

        Switch inScatter_;
        dictionary coeffsDict_;
        label subAngleNum_;

        //- Particle radius [m]
        scalar radius_;

        //- Complex particle refractive index
        std::complex<scalar> mParticle_;

        //- Real refractive index of the host medium
        scalar mMedium_;

        //- Per-band central wavelengths in vacuum [nm]
        List<scalar> wavelengths_;

        //- Pixel-integrated, row-normalised phase function table.
        //  Indexed phaseFunction_[j + i*nAngle + iband*nAngle*nAngle].
        DynamicList<scalar, 1, 0> phaseFunction_;

public:

    //- Runtime type information
    TypeName("mieModel");


    // Constructors

        //- Construct from components
        mieModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~mieModel();


    // Member Functions

        virtual bool inScatter() const
        {
            return inScatter_;
        }

        Foam::scalar correct
        (
            const label rayI,
            const label rayJ,
            const label iBand
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
