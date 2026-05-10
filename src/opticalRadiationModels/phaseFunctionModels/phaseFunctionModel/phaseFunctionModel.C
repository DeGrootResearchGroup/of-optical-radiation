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

#include "phaseFunctionModel.H"
#include "DOM.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(phaseFunctionModel, 0);
        defineRunTimeSelectionTable(phaseFunctionModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::phaseFunctionModel::phaseFunctionModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    dom_(dom),
    dict_(dict),
    nDim_(nDim),
    nBand_(0),
    nAngle_(0),
    inScatter_(false),
    subAngleNum_(1),
    phaseFunction_()
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::optical::phaseFunctionModel::~phaseFunctionModel()
{}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::optical::phaseFunctionModel::buildPhaseTable()
{
    // Angular grid:
    //   nPhi azimuthal divisions per pi (so 2*nPhi cover [0, 2 pi)).
    //   nTheta polar divisions over [0, pi].
    //   deltaPhi   = pi / (2 nPhi)
    //   deltaTheta = pi / nTheta
    // Each Omega_j is a (deltaPhi by deltaTheta) bin centred at the j-th
    // ray's (phi, theta).
    const label nPhi   = dom_.nPhi();
    const label nTheta = dom_.nTheta();
    const scalar deltaPhi   = pi/(2.0*nPhi);
    const scalar deltaTheta = pi/nTheta;

    // Sub-pixel sizes inside one Omega_j.
    const scalar dp = deltaPhi/subAngleNum_;
    const scalar dt = deltaTheta/subAngleNum_;

    phaseFunction_.setSize(nAngle_*nAngle_*nBand_);
    forAll(phaseFunction_, k)
    {
        phaseFunction_[k] = 0.0;
    }

    for (label iBand = 0; iBand < nBand_; ++iBand)
    {
        for (label i = 0; i < nAngle_; ++i)
        {
            const label rayI = i + iBand*nAngle_;
            const vector dI = dom_.IRay(rayI).d();

            // Pass 1: pixel-sample Phi over each Omega_j to fill Psi_ij,
            // and accumulate pfSum_i = sum_j Psi_ij in parallel.
            scalar pfSum = 0;
            for (label j = 0; j < nAngle_; ++j)
            {
                const label rayJ = j + iBand*nAngle_;
                const label idx = j + i*nAngle_ + iBand*nAngle_*nAngle_;

                const scalar phiJ   = dom_.IRay(rayJ).phi();
                const scalar thetaJ = dom_.IRay(rayJ).theta();

                if (nDim_ == 3)
                {
                    // 3-D: subAngleNum x subAngleNum sub-pixels in
                    // (phi, theta), each with solid angle
                    //   nOmega = 2 sin(theta_pixel) sin(dt/2) dp
                    // (the sin(theta) Jacobian comes from dOmega = sin
                    // theta dtheta dphi).
                    for (label m = 0; m < subAngleNum_; ++m)
                    {
                        for (label n = 0; n < subAngleNum_; ++n)
                        {
                            const scalar nP =
                                phiJ + (2.0*m - subAngleNum_ + 1.0)*dp/2.0;
                            const scalar nT =
                                thetaJ
                              + (2.0*n - subAngleNum_ + 1.0)*dt/2.0;
                            const scalar nOmega =
                                2.0*sin(nT)*sin(dt/2.0)*dp;
                            const vector nD
                            (
                                sin(nT)*cos(nP),
                                sin(nT)*sin(nP),
                                cos(nT)
                            );
                            const scalar cosV = dI & nD;
                            phaseFunction_[idx] +=
                                phaseShape(cosV, iBand)*nOmega;
                        }
                    }
                }
                else // nDim_ == 2: rays only in the xy-plane
                {
                    // 2-D: subAngleNum sub-pixels in phi only; the polar
                    // direction is handled by the 2-D-as-3-D scheme that
                    // sets theta = pi/2 and deltaTheta = pi for the rays.
                    // The pixel solid-angle weight 2*dp matches OF's
                    // 2-D ray omega = 2*deltaPhi.
                    for (label m = 0; m < subAngleNum_; ++m)
                    {
                        const scalar nP =
                            phiJ + (2.0*m - subAngleNum_ + 1.0)*dp/2.0;
                        const scalar nOmega = 2.0*dp;
                        const vector nD(cos(nP), sin(nP), 0);
                        const scalar cosV = dI & nD;
                        phaseFunction_[idx] +=
                            phaseShape(cosV, iBand)*nOmega;
                    }
                }

                pfSum += phaseFunction_[idx];
            }

            // Pass 2: row-normalise so sum_j table[i, j, iBand] = 1.
            // pfSum_i is the discrete approximation of integral over
            // 4 pi of phaseShape dOmega; for a phaseShape that already
            // integrates to 1 (HG with the (1-g^2)/(4 pi) prefactor) the
            // normalisation is a small correction for pixel
            // discretisation error, while for an unscaled shape (e.g.
            // raw BHMIE intensity) it provides the absolute scaling.
            for (label j = 0; j < nAngle_; ++j)
            {
                const label idx = j + i*nAngle_ + iBand*nAngle_*nAngle_;
                phaseFunction_[idx] /= pfSum;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::optical::phaseFunctionModel::correct
(
    const label rayI,
    const label rayJ,
    const label iBand
) const
{
    if (phaseFunction_.empty())
    {
        return 0.0;
    }
    const label i = rayI - iBand*nAngle_;
    const label j = rayJ - iBand*nAngle_;
    return phaseFunction_[j + i*nAngle_ + iBand*nAngle_*nAngle_];
}


// ************************************************************************* //
