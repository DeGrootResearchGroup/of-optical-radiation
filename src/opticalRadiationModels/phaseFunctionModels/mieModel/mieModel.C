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
#include "mieKernel.H"
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
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    radius_(0),
    mParticle_(0, 0),
    mMedium_(1),
    wavelengths_(),
    phaseFunction_()
{
    coeffsDict_.lookup("inScatter") >> inScatter_;

    if (!inScatter_)
    {
        return;
    }

    nBand_  = dom_.nBand();
    nAngle_ = dom_.nAngle();
    coeffsDict_.lookup("subAngleNum") >> subAngleNum_;

    radius_  = readScalar(coeffsDict_.lookup("radius"));
    mMedium_ = coeffsDict_.lookupOrDefault<scalar>("mMedium", 1.0);
    {
        Vector2D<scalar> mp;
        coeffsDict_.lookup("mParticle") >> mp;
        mParticle_ = std::complex<scalar>(mp.x(), mp.y());
    }
    coeffsDict_.lookup("wavelengths") >> wavelengths_;

    if (wavelengths_.size() != nBand_)
    {
        FatalErrorIn
        (
            "mieModel::mieModel(const DOM&, const dictionary&, const label&)"
        )   << "wavelengths list has " << wavelengths_.size()
            << " entries but DOM nBand = " << nBand_
            << exit(FatalError);
    }
    if (radius_ <= 0)
    {
        FatalErrorIn
        (
            "mieModel::mieModel(const DOM&, const dictionary&, const label&)"
        )   << "radius must be positive (got " << radius_ << " m)"
            << exit(FatalError);
    }
    if (mMedium_ <= 0)
    {
        FatalErrorIn
        (
            "mieModel::mieModel(const DOM&, const dictionary&, const label&)"
        )   << "mMedium must be positive (got " << mMedium_ << ")"
            << exit(FatalError);
    }

    const std::complex<scalar> mRel(mParticle_/mMedium_);

    const label nPhi   = dom_.nPhi();
    const label nTheta = dom_.nTheta();
    const scalar deltaPhi   = pi/(2.0*nPhi);
    const scalar deltaTheta = pi/nTheta;

    phaseFunction_.setSize(nAngle_*nAngle_*nBand_);
    forAll(phaseFunction_, i)
    {
        phaseFunction_[i] = 0.0;
    }

    // Per band, build the Mie kernel once and reuse it for every (i,j)
    // and every sub-pixel evaluation. Kernel construction is O(N^2) where
    // N = ceil(x + 4 x^(1/3) + 2); a phaseIntensity() call is O(N).
    for (label iBand = 0; iBand < nBand_; ++iBand)
    {
        const scalar lamM = wavelengths_[iBand]*1e-9;
        const scalar x = 2.0*pi*radius_*mMedium_/lamM;

        mieKernel mie(x, mRel);

        Info<< "    mie phase function band " << iBand
            << ": lambda = " << wavelengths_[iBand] << " nm, "
            << "x = " << x
            << ", g = " << mie.g()
            << ", N_terms = " << mie.N() << endl;

        if (nDim == 3)
        {
            const scalar dp = deltaPhi/subAngleNum_;
            const scalar dt = deltaTheta/subAngleNum_;

            for (label i = 0; i < nAngle_; ++i)
            {
                const label rayI = i + iBand*nAngle_;
                scalar pfSum = 0;
                for (label j = 0; j < nAngle_; ++j)
                {
                    const label rayJ = j + iBand*nAngle_;
                    const label idx =
                        j + i*nAngle_ + iBand*nAngle_*nAngle_;
                    for (label m = 0; m < subAngleNum_; ++m)
                    {
                        for (label n = 0; n < subAngleNum_; ++n)
                        {
                            const scalar nP =
                                (2.0*m - subAngleNum_ + 1.0)*dp/2.0
                              + dom.IRay(rayJ).phi();
                            const scalar nT =
                                (2.0*n - subAngleNum_ + 1.0)*dt/2.0
                              + dom.IRay(rayJ).theta();
                            const scalar nOmega =
                                2*sin(nT)*sin(dt/2)*dp;
                            const vector nD = vector
                            (
                                sin(nT)*cos(nP),
                                sin(nT)*sin(nP),
                                cos(nT)
                            );
                            const scalar cosV = dom.IRay(rayI).d() & nD;
                            phaseFunction_[idx] +=
                                mie.phaseIntensity(cosV)*nOmega;
                        }
                    }
                    pfSum += phaseFunction_[idx];
                    phaseFunction_[idx] /= dom.IRay(rayI).omega();
                }

                for (label j = 0; j < nAngle_; ++j)
                {
                    const label idx =
                        j + i*nAngle_ + iBand*nAngle_*nAngle_;
                    phaseFunction_[idx] /= pfSum;
                }
            }
        }

        if (nDim == 2)
        {
            const scalar dp = deltaPhi/subAngleNum_;

            for (label i = 0; i < nAngle_; ++i)
            {
                const label rayI = i + iBand*nAngle_;
                scalar pfSum = 0;
                for (label j = 0; j < nAngle_; ++j)
                {
                    const label rayJ = j + iBand*nAngle_;
                    const label idx =
                        j + i*nAngle_ + iBand*nAngle_*nAngle_;
                    for (label m = 0; m < subAngleNum_; ++m)
                    {
                        const scalar nP =
                            (2.0*m - subAngleNum_ + 1.0)*dp/2.0
                          + dom.IRay(rayJ).phi();
                        const scalar nOmega = 2*dp;
                        const vector nD = vector(cos(nP), sin(nP), 0);
                        const scalar cosV = dom.IRay(rayI).d() & nD;
                        phaseFunction_[idx] +=
                            mie.phaseIntensity(cosV)*nOmega;
                    }
                    pfSum += phaseFunction_[idx];
                    phaseFunction_[idx] /= dom.IRay(rayI).omega();
                }

                for (label j = 0; j < nAngle_; ++j)
                {
                    const label idx =
                        j + i*nAngle_ + iBand*nAngle_*nAngle_;
                    phaseFunction_[idx] /= pfSum;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::mieModel::~mieModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::optical::mieModel::correct
(
    const label rayI,
    const label rayJ,
    const label iBand
) const
{
    // rayI / rayJ are flat indices that already include the band offset
    // (rayI = iAngleI + iBand*nAngle_); the table is laid out by per-band
    // angle indices, so subtract the offset before forming the index.
    const label i = rayI - iBand*nAngle_;
    const label j = rayJ - iBand*nAngle_;
    return phaseFunction_[j + i*nAngle_ + iBand*nAngle_*nAngle_];
}


// ************************************************************************* //
