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

\*---------------------------------------------------------------------------*/

#include "mieModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"
#include "Vector2D.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(mieModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            mieModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::mieModel::mieModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom, dict, nDim),
    radius_(0),
    mParticle_(0, 0),
    mMedium_(1),
    wavelengths_(),
    kernels_()
{
    const dictionary& coeffs = dict.subDict(typeName + "Coeffs");
    coeffs.lookup("inScatter") >> inScatter_;

    if (!inScatter_)
    {
        return;
    }

    nBand_  = dom_.nBand();
    nAngle_ = dom_.nAngle();
    coeffs.lookup("subAngleNum") >> subAngleNum_;

    radius_  = readScalar(coeffs.lookup("radius"));
    mMedium_ = coeffs.lookupOrDefault<scalar>("mMedium", 1.0);
    {
        Vector2D<scalar> mp;
        coeffs.lookup("mParticle") >> mp;
        mParticle_ = std::complex<scalar>(mp.x(), mp.y());
    }
    coeffs.lookup("wavelengths") >> wavelengths_;

    if (wavelengths_.size() != nBand_)
    {
        FatalErrorInFunction
            << "wavelengths list has " << wavelengths_.size()
            << " entries but DOM nBand = " << nBand_
            << exit(FatalError);
    }
    if (radius_ <= 0)
    {
        FatalErrorInFunction
            << "radius must be positive (got " << radius_ << " m)"
            << exit(FatalError);
    }
    if (mMedium_ <= 0)
    {
        FatalErrorInFunction
            << "mMedium must be positive (got " << mMedium_ << ")"
            << exit(FatalError);
    }

    // Build one kernel per band. Kernel construction is O(N^2) where
    // N = ceil(x + 4 x^(1/3) + 2); a phaseIntensity() call is O(N).
    const std::complex<scalar> mRel(mParticle_/mMedium_);
    kernels_.setSize(nBand_);
    for (label iBand = 0; iBand < nBand_; ++iBand)
    {
        const scalar lamM = wavelengths_[iBand]*1e-9;
        const scalar x = 2.0*pi*radius_*mMedium_/lamM;

        kernels_.set(iBand, new mieKernel(x, mRel));

        Info<< "    mie phase function band " << iBand
            << ": lambda = " << wavelengths_[iBand] << " nm, "
            << "x = " << x
            << ", g = " << kernels_[iBand].g()
            << ", N_terms = " << kernels_[iBand].N() << endl;
    }

    buildPhaseTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::mieModel::~mieModel()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::optical::mieModel::phaseShape
(
    const scalar cosV,
    const label iBand
) const
{
    return kernels_[iBand].phaseIntensity(cosV);
}


// ************************************************************************* //
