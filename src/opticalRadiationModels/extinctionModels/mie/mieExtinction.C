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

#include "mieExtinction.H"
#include "mieKernel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(mieExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            mieExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::mieExtinction::mieExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    nBands_(readLabel(coeffsDict_.lookup("nBands"))),
    wavelengths_(),
    radius_(readScalar(coeffsDict_.lookup("radius"))),
    mParticle_(0, 0),
    mMedium_(coeffsDict_.lookupOrDefault<scalar>("mMedium", 1.0)),
    numberDensityFieldName_
    (
        coeffsDict_.lookupOrDefault<word>("numberDensityField", "n")
    ),
    QscaPerBand_(nBands_, 0.0),
    QabsPerBand_(nBands_, 0.0),
    sigmaScaPrefactorPerBand_(nBands_, 0.0),
    sigmaAbsPrefactorPerBand_(nBands_, 0.0)
{
    init(nBands_);

    coeffsDict_.lookup("wavelengths") >> wavelengths_;

    if (wavelengths_.size() != nBands_)
    {
        FatalErrorIn
        (
            "mieExtinction::mieExtinction"
            "(const dictionary&, const fvMesh&)"
        )   << "wavelengths list has " << wavelengths_.size()
            << " entries but nBands = " << nBands_
            << exit(FatalError);
    }

    {
        Vector2D<scalar> mp;
        coeffsDict_.lookup("mParticle") >> mp;
        mParticle_ = std::complex<scalar>(mp.x(), mp.y());
    }

    if (radius_ <= 0)
    {
        FatalErrorIn
        (
            "mieExtinction::mieExtinction"
            "(const dictionary&, const fvMesh&)"
        )   << "radius must be positive (got " << radius_ << " m)"
            << exit(FatalError);
    }
    if (mMedium_ <= 0)
    {
        FatalErrorIn
        (
            "mieExtinction::mieExtinction"
            "(const dictionary&, const fvMesh&)"
        )   << "mMedium must be positive (got " << mMedium_ << ")"
            << exit(FatalError);
    }

    const scalar pi = constant::mathematical::pi;

    // Geometric cross-section pi r^2 [m^2] -- folded into the per-band
    // prefactor so correct() reduces to a single multiplication by N(x).
    const scalar piR2 = pi*radius_*radius_;

    // Relative refractive index (BHMIE convention: m_particle / m_medium,
    // with imag(m_particle) >= 0 for an absorbing material).
    const std::complex<scalar> mRel(mParticle_/mMedium_);

    forAll(wavelengths_, b)
    {
        const scalar lamM = wavelengths_[b]*1e-9;     // nm -> m
        const scalar x = 2.0*pi*radius_*mMedium_/lamM;

        mieKernel mie(x, mRel);

        QscaPerBand_[b] = mie.Qsca();
        // Q_abs = Q_ext - Q_sca can pick up a small negative roundoff
        // for purely scattering particles; clip at zero so a non-
        // physical negative absorption coefficient never enters the
        // RTE. Floating-point noise is the only thing affected.
        QabsPerBand_[b] = max(scalar(0), mie.Qabs());

        sigmaScaPrefactorPerBand_[b] = piR2*QscaPerBand_[b];
        sigmaAbsPrefactorPerBand_[b] = piR2*QabsPerBand_[b];

        Info<< "    mie band " << b
            << ": lambda = " << wavelengths_[b] << " nm, "
            << "x = " << x
            << ", Q_sca = " << mie.Qsca()
            << ", Q_abs = " << mie.Qabs()
            << ", g = " << mie.g()
            << ", N_terms = " << mie.N() << endl;
    }

    // Load the number-density field from disk if no host solver has
    // already registered it.
    if (!mesh.foundObject<volScalarField>(numberDensityFieldName_))
    {
        ownedNumberDensityField_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    numberDensityFieldName_,
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::mieExtinction::~mieExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::mieExtinction::correct()
{
    const volScalarField& N =
        mesh_.lookupObject<volScalarField>(numberDensityFieldName_);

    forAll(ALambda_, iBand)
    {
        const dimensionedScalar sigmaA
        (
            "sigmaA",
            dimLength*dimLength,
            sigmaAbsPrefactorPerBand_[iBand]
        );
        const dimensionedScalar sigmaS
        (
            "sigmaS",
            dimLength*dimLength,
            sigmaScaPrefactorPerBand_[iBand]
        );

        ALambda_[iBand] = sigmaA*N;
        SLambda_[iBand] = sigmaS*N;
    }
}


// ************************************************************************* //
