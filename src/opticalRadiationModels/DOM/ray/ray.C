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

#include "ray.H"
#include "fvm.H"
#include "DOM.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word
Foam::optical::ray::namePrefix("I");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::ray::ray
(
    const DOM& dom,
    const fvMesh& mesh,
    const label iBand,
    const label iAngle,
    const scalar theta,
    const scalar phi,
    const scalar deltaTheta,
    const scalar deltaPhi
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

    dAve_ = dom_.intDirOmega(theta_, phi_);

    autoPtr<volScalarField> IDefaultPtr;

    IOobject IHeader
    (
        namePrefix + "_" + name(iBand) + "_" + name(iAngle),
        mesh_.time().name(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    // check if field exists and can be read
    if (IHeader.headerOk())
    {
        I_.reset(new volScalarField(IHeader, mesh_));
    }
    else
    {
        // Demand driven load the IDefault field
        if (IDefaultPtr.empty())
        {
            IDefaultPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "I",
                        mesh_.time().name(),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::ray::~ray()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::optical::ray::correct(const volScalarField& ds)
{
    scalar maxResidual = -GREAT;
    scalar eqnResidual;

    const volScalarField& A = dom_.A(iBand_);
    const volScalarField& S = dom_.S(iBand_);
    const volScalarField K = A + S;

    surfaceScalarField Ji0(vector(0,0,0) & mesh_.Sf());
    surfaceScalarField Ji1(vector(0,0,0) & mesh_.Sf());
    computeFluxCoeffs_(Ji0, Ji1);

    fvScalarMatrix IiEq
    (
          fvm::div(Ji0, I_(), "div(Ji,Ii_h)")
        + fvm::div(Ji1, I_(), "div(Ji,Ii_h)")
        + fvm::Sp(K*omega_, I_())
       == S*ds*omega_
    );

    IiEq.relax();
    eqnResidual = solve(IiEq, word("Ii")).initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    return maxResidual;
}


void Foam::optical::ray::computeFluxCoeffs_
(
    surfaceScalarField& Ji0,
    surfaceScalarField& Ji1
) const
{
    const label npTheta = dom_.nPixelTheta();
    const label npPhi   = dom_.nPixelPhi();
    const scalar deltaTheta  = dom_.deltaTheta();
    const scalar deltaPhi    = dom_.deltaPhi();
    const scalar pixelDTheta = deltaTheta/npTheta;
    const scalar pixelDPhi   = deltaPhi  /npPhi;

    const surfaceVectorField& Sf = mesh_.Sf();
    const vectorField& Sf_int = Sf.primitiveField();
    scalarField& Ji0_int = Ji0.primitiveFieldRef();
    scalarField& Ji1_int = Ji1.primitiveFieldRef();

    const surfaceVectorField::Boundary& Sf_bf = Sf.boundaryField();
    surfaceScalarField::Boundary& Ji0_bf = Ji0.boundaryFieldRef();
    surfaceScalarField::Boundary& Ji1_bf = Ji1.boundaryFieldRef();

    for (label i = 0; i < npTheta; i++)
    {
        for (label j = 0; j < npPhi; j++)
        {
            const scalar pixelTheta =
                theta_ - 0.5*deltaTheta + (i + 0.5)*pixelDTheta;
            const scalar pixelPhi =
                phi_   - 0.5*deltaPhi   + (j + 0.5)*pixelDPhi;

            const vector pixelDir =
                dom_.anglesToDir(pixelTheta, pixelPhi);
            const vector pixelFlux =
                dom_.intDirOmega
                (
                    pixelTheta, pixelPhi, pixelDTheta, pixelDPhi
                );

            forAll(Sf_int, fi)
            {
                const scalar dpd = pixelDir & Sf_int[fi];
                if (dpd > 0)
                {
                    Ji0_int[fi] += pixelFlux & Sf_int[fi];
                }
                else if (dpd < 0)
                {
                    Ji1_int[fi] += pixelFlux & Sf_int[fi];
                }
            }

            forAll(Sf_bf, patchi)
            {
                const fvsPatchField<vector>& Sfp  = Sf_bf[patchi];
                fvsPatchField<scalar>&       Ji0p = Ji0_bf[patchi];
                fvsPatchField<scalar>&       Ji1p = Ji1_bf[patchi];

                forAll(Sfp, fi)
                {
                    const scalar dpd = pixelDir & Sfp[fi];
                    if (dpd > 0)
                    {
                        Ji0p[fi] += pixelFlux & Sfp[fi];
                    }
                    else if (dpd < 0)
                    {
                        Ji1p[fi] += pixelFlux & Sfp[fi];
                    }
                }
            }
        }
    }
}


void Foam::optical::ray::updateBoundary()
{
    I_->correctBoundaryConditions();
}

// ************************************************************************* //
