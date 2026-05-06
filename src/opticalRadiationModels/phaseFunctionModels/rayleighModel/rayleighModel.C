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

#include "rayleighModel.H"
#include "addToRunTimeSelectionTable.H"
#include "DOM.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(rayleighModel, 0);
        addToRunTimeSelectionTable
        (
            phaseFunctionModel,
            rayleighModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::rayleighModel::rayleighModel
(
    const DOM& dom,
    const dictionary& dict,
    const label& nDim
)
:
    phaseFunctionModel(dom, dict, nDim),
    coeffsDict_(dict.subDict(typeName + "Coeffs"))
{
    coeffsDict_.lookup("inScatter") >> inScatter_;

    if (inScatter_)
    {
        nBand_ = dom_.nBand();
        nAngle_ = dom_.nAngle();
        coeffsDict_.lookup("subAngleNum") >> subAngleNum_;

        const label nPhi   = dom_.nPhi();
        const label nTheta = dom_.nTheta();
        const scalar deltaPhi   = pi/(2.0*nPhi);
        const scalar deltaTheta = pi/nTheta;

        phaseFunction_.setSize(nAngle_*nAngle_*nBand_);
        for (label i = 0; i < nAngle_*nAngle_*nBand_; i++)
        {
            phaseFunction_[i] = 0.0;
        }

        if (nDim == 3)
        {
            const scalar dp = deltaPhi/subAngleNum_;
            const scalar dt = deltaTheta/subAngleNum_;

            for (label iband = 0; iband < nBand_; iband++)
            {
                for (label i = 0; i < nAngle_; i++)
                {
                    const label rayI = i + iband*nAngle_;
                    scalar pfSum = 0;
                    for (label j = 0; j < nAngle_; j++)
                    {
                        const label rayJ = j + iband*nAngle_;
                        const label idx =
                            j + i*nAngle_ + iband*nAngle_*nAngle_;
                        for (label m = 0; m < subAngleNum_; m++)
                        {
                            for (label n = 0; n < subAngleNum_; n++)
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
                                const scalar cosV =
                                    dom.IRay(rayI).d() & nD;

                                phaseFunction_[idx] +=
                                    rayleigh3d(cosV)*nOmega;
                            }
                        }
                        pfSum += phaseFunction_[idx];
                        phaseFunction_[idx] /= dom.IRay(rayI).omega();
                    }

                    for (label j = 0; j < nAngle_; j++)
                    {
                        const label idx =
                            j + i*nAngle_ + iband*nAngle_*nAngle_;
                        phaseFunction_[idx] /= pfSum;
                    }
                }
            }
        }

        if (nDim == 2)
        {
            const scalar dp = deltaPhi/subAngleNum_;

            for (label iband = 0; iband < nBand_; iband++)
            {
                for (label i = 0; i < nAngle_; i++)
                {
                    const label rayI = i + iband*nAngle_;
                    scalar pfSum = 0;
                    for (label j = 0; j < nAngle_; j++)
                    {
                        const label rayJ = j + iband*nAngle_;
                        const label idx =
                            j + i*nAngle_ + iband*nAngle_*nAngle_;
                        for (label m = 0; m < subAngleNum_; m++)
                        {
                            const scalar nP =
                                (2.0*m - subAngleNum_ + 1.0)*dp/2.0
                              + dom.IRay(rayJ).phi();
                            const scalar nOmega = 2*dp;
                            const vector nD = vector
                            (
                                cos(nP), sin(nP), 0
                            );
                            const scalar cosV =
                                dom.IRay(rayI).d() & nD;

                            phaseFunction_(idx) +=
                                rayleigh2d(cosV)*nOmega;
                        }
                        pfSum += phaseFunction_(idx);
                        phaseFunction_[idx] /= dom.IRay(rayI).omega();
                    }

                    for (label j = 0; j < nAngle_; j++)
                    {
                        const label idx =
                            j + i*nAngle_ + iband*nAngle_*nAngle_;
                        phaseFunction_[idx] /= pfSum;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::rayleighModel::~rayleighModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::optical::rayleighModel::correct
(
    const label rayI,
    const label rayJ,
    const label iBand
) const
{
    // Same per-band index reconstruction as HenyeyGreensteinModel: rayI/rayJ
    // are flat indices that include the band offset; the table is laid out
    // by per-band angle indices.
    const label i = rayI - iBand*nAngle_;
    const label j = rayJ - iBand*nAngle_;
    return phaseFunction_[j + i*nAngle_ + iBand*nAngle_*nAngle_];
}


Foam::scalar Foam::optical::rayleighModel::rayleigh3d
(
    const scalar cosV
) const
{
    return (3.0/(16.0*pi))*(1.0 + cosV*cosV);
}


Foam::scalar Foam::optical::rayleighModel::rayleigh2d
(
    const scalar cosV
) const
{
    return (1.0 + cosV*cosV)/(3.0*pi);
}


// ************************************************************************* //
