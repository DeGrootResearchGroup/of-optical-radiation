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

#include "photoBioIntensityRay.H"
#include "fvm.H"
#include "photoBioDOM.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::photoBio::photoBioIntensityRay::rayId(0);

const Foam::word
Foam::photoBio::photoBioIntensityRay::intensityPrefix("I");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::photoBioIntensityRay::photoBioIntensityRay
(
    const photoBioDOM& dom,
    const fvMesh& mesh,
    const label iBand,
    const label iAngle,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta
)
:
    dom_(dom),
    mesh_(mesh),
    d_(vector::zero),
    dAve_(vector::zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    iBand_(iBand),
    iAngle_(iAngle)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;

    d_ = vector(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
    dAve_ = vector
    (
       cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi
        *Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );

    autoPtr<volScalarField> IDefaultPtr;

    IOobject IHeader
    (
        intensityPrefix + "_" + name(iBand) + "_" + name(iAngle),
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    // check if field exists and can be read
    if (IHeader.typeHeaderOk<volScalarField>())
    {
        I_.reset(new volScalarField(IHeader, mesh_));
    }
    else
    {
        // Demand driven load the IDefault field
        if (!IDefaultPtr)
        {
            IDefaultPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "I",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_
                )
            );
        }

        // Construct header with NO_READ so the field copies IDefault values
        IOobject noReadHeader
        (
            IHeader.name(),
            IHeader.instance(),
            IHeader.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        I_.reset(new volScalarField(noReadHeader, IDefaultPtr()));
    }
    rayId++;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::photoBio::photoBioIntensityRay::~photoBioIntensityRay()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::photoBio::photoBioIntensityRay::correct()
{
    scalar maxResidual = -GREAT;
    scalar eqnResidual;

    const volScalarField& A = dom_.A(iBand_);
    const volScalarField& S = dom_.S(iBand_);
    const volScalarField K = A + S;

    const volScalarField& ds = dom_.diffusionScatter();
    const surfaceScalarField Ji = dAve_ & mesh_.Sf();

    fvScalarMatrix IiEq
    (
          fvm::div(Ji, I_(), "div(Ji,Ii_h)")
        + fvm::Sp(K*omega_, I_())
       == S*ds*omega_
    );

    IiEq.relax();
    eqnResidual = solve(IiEq, mesh_.solver("Ii")).initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    return maxResidual;
}


void Foam::photoBio::photoBioIntensityRay::updateBoundary()
{
    I_->correctBoundaryConditions();
}

// ************************************************************************* //
