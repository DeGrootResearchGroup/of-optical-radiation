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

\*---------------------------------------------------------------------------*/

inline const Foam::optical::ray&
Foam::optical::DOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline Foam::label Foam::optical::DOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::optical::DOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::optical::DOM::nAngle() const
{
    return nAngle_;
}


inline Foam::label Foam::optical::DOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::optical::DOM::nBand() const
{
    return nBand_;
}


inline Foam::label Foam::optical::DOM::nPixelPhi() const
{
    return nPixelPhi_;
}


inline Foam::label Foam::optical::DOM::nPixelTheta() const
{
    return nPixelTheta_;
}


inline Foam::scalar Foam::optical::DOM::deltaTheta() const
{
    return deltaTheta_;
}


inline Foam::scalar Foam::optical::DOM::deltaPhi() const
{
    return deltaPhi_;
}


inline const Foam::volScalarField& Foam::optical::DOM::A(const label iLambda) const
{
    return extinction_->A(iLambda);
}


inline const Foam::volScalarField& Foam::optical::DOM::S(const label iLambda) const
{
    return extinction_->S(iLambda);
}


inline const Foam::volScalarField& Foam::optical::DOM::G() const
{
    return G_;
}


// ************************************************************************* //
