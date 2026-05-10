/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
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

#include "DOM.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "phaseFunctionModel.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(DOM, 0);
        addToRunTimeSelectionTable
        (
            radiationModel,
            DOM,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::DOM::DOM(const volScalarField& I)
:
    radiationModel(typeName, I),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        // [W/m^2] = kg/s^3
        dimensionedScalar("G", dimPower/dimArea, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nAngle_(0),
    nRay_(0),
    nBand_(coeffs_.lookupOrDefault<label>("nBand", 1)),
    nPixelPhi_(coeffs_.lookupOrDefault<label>("nPixelPhi", 1)),
    nPixelTheta_(coeffs_.lookupOrDefault<label>("nPixelTheta", 1)),
    GLambda_(nBand_),
    IRay_(0),
    ISnapshot_(0),
    convergence_(readScalar(coeffs_.lookup("convergence"))),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50))
{
    Info<< "Creating DOM model with " << nBand_ << " bands" << endl;

    // Check that dimension of mesh is compatible with settings
    checkDim_();

    // Set the number of angles (=2*nPhi*nTheta)
    nAngle_ = 2*nPhi_*nTheta_;

    // Set the number of rays, and allocate pointer list
    nRay_ = nAngle_*nBand_;
    IRay_.setSize(nRay_);

    // Set deltaTheta (=pi/nTheta for 3D; =pi for 2D)
    deltaTheta_ = mesh_.nSolutionD() == 3 ? pi/nTheta_ : pi;

    // Set deltaPhi
    deltaPhi_ = pi/nPhi_;

    // Set up all of the rays
    label i = 0;
    for (label iBand = 0; iBand < nBand_; iBand++)
    {
        for (label iTheta = 0; iTheta < nTheta_; iTheta++)
        {
            for (label iPhi = 0; iPhi < nAngle_/nTheta_; iPhi++)
            {
                label iAngle = iPhi + 2*iTheta*nPhi_;
                scalar theta = (iTheta + 0.5)*deltaTheta_;
                scalar phi = (iPhi + 0.5)*deltaPhi_;
                setRay_(i, iBand, iAngle, theta, phi);
                i++;
            }
        }
    }

    // Verify the angular discretisation covers the full sphere; sumOmega
    // is band-independent so we only check the first band's rays.
    scalar sumOmega = 0.0;
    for (label rayI = 0; rayI < nAngle_; rayI++)
    {
        sumOmega += IRay_[rayI].omega();
    }
    Info << "Sum of solid angles: " << sumOmega
        << " (expected " << 4.0*pi << ")" << endl;
    const scalar omegaTol = 1e-6 * 4.0*pi;
    if (mag(sumOmega - 4.0*pi) > omegaTol)
    {
        WarningInFunction
            << "Sum of solid angles deviates from 4*pi by "
            << mag(sumOmega - 4.0*pi)
            << ", which exceeds the tolerance of " << omegaTol << "."
            << " The angular discretisation may be malformed." << endl;
    }

    Info<< "DOM : Allocated " << IRay_.size() << " rays" << endl;

    // Allocate the per-ray radiance snapshot (one volScalarField per ray,
    // matching IRay_[i].I() in size and dimensions). Used by calculate()
    // to freeze I_j values at the start of each outer iteration so the
    // in-scatter source uses a Jacobi update.
    ISnapshot_.setSize(nRay_);
    forAll(IRay_, rayI)
    {
        ISnapshot_.set
        (
            rayI,
            new volScalarField
            (
                IOobject
                (
                    "ISnapshot_" + Foam::name(rayI),
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                IRay_[rayI].I()
            )
        );
    }

    forAll(GLambda_, iBand)
    {
        GLambda_.set
        (
            iBand,
            new volScalarField
            (
                IOobject
                (
                    "GLambda_" + Foam::name(iBand) ,
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                G_
            )
        );
    }

    phaseFunctionModel_ = phaseFunctionModel::New(*this,coeffs_, mesh_.nSolutionD());

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::DOM::~DOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::optical::DOM::read()
{
    if (radiationModel::read())
    {
        // Solution parameters. Ray geometry (nTheta, nPhi, nBand,
        // pixelation) is held constant -- changing it would require
        // tearing down and rebuilding the per-ray storage, which is
        // out of scope for a runtime-modifiable read.
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        // Refresh the phase function model the same way the base
        // class refreshes the extinction model: concrete phase-
        // function models bake their dictionary into per-band
        // tables at construction and don't expose a read(), so a
        // factory rebuild is the way to pick up coefficient edits
        // (e.g. asymmetry parameter g for HG / Schlick).
        phaseFunctionModel_.clear();
        phaseFunctionModel_ =
            phaseFunctionModel::New(*this, coeffs_, mesh_.nSolutionD());

        return true;
    }
    else
    {
        return false;
    }
}

void Foam::optical::DOM::calculate()
{
    // Correct the extinction model.
    extinction_->correct();

    // In-scatter source S_in,i scratch field, allocated once per
    // calculate(). Re-zeroed and re-filled per (iBand, rayI) pair
    // inside the inner loop. [W/m^2 = kg/s^3] -- same dimensionSet
    // group as I, modulo the implicit /sr that DOM tracks
    // dimensionally as dimensionless throughout.
    volScalarField ds
    (
        IOobject
        (
            "diffusionScatter",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("diffusionScatter", dimPower/dimArea, 0.0)
    );

    const bool doInScatter = phaseFunctionModel_->inScatter();
    scalar maxResidual = 0.0;
    label radIter = 0;

    do
    {
        radIter++;
        maxResidual = 0.0;

        // Snapshot the current per-ray radiance fields for use as the
        // in-scatter source. The inner loop below updates IRay_[i].I()
        // sequentially (Gauss-Seidel over rays), and using the partially
        // updated values to build ds for ray i can drive an oscillation
        // in the outer iteration on strongly-coupled cases. Computing ds
        // from a snapshot frozen at the start of the outer iteration
        // (Jacobi) symmetrises the coupling and stabilises convergence.
        if (doInScatter)
        {
            forAll(IRay_, rayI)
            {
                ISnapshot_[rayI] == IRay_[rayI].I();
            }
        }

        forAll(IRay_, rayI)
        {
            const label iBand = IRay_[rayI].iBand();

            if (debug)
            {
                Info<< "opticalRadiation solver:"
                    << "    iter: "   << radIter
                    << "    iBand: "  << iBand
                    << "    iAngle: " << IRay_[rayI].iAngle() << endl;
            }

            if (doInScatter)
            {
                // Fused in-scatter accumulator. The table value
                // table[i, j, iBand] from buildPhaseTable is row-
                // normalised so sum_j table[i, j, iBand] = 1, which
                // means it absorbs both the per-bin solid angle
                // omega_j and the 1/(4 pi) prefactor of the in-scatter
                // integral. The discrete sum here approximates
                //   (1/(4 pi)) integral over 4 pi of
                //                Phi(s_i . s') I(s') dOmega'
                // with no extra solid-angle weight.
                //
                // The whole row table[i, *, iBand] is fetched once via
                // phaseRow() (one virtual call per (rayI, iBand) pair),
                // and the cell- and boundary-loops below fold the
                // per-pair multiply-and-accumulate into a single pass
                // -- no tmp<volScalarField> per ray pair.
                const scalar* row =
                    phaseFunctionModel_->phaseRow(rayI, iBand);

                scalarField& dsCell = ds.primitiveFieldRef();
                dsCell = 0.0;

                volScalarField::Boundary& dsBf = ds.boundaryFieldRef();
                forAll(dsBf, patchi)
                {
                    dsBf[patchi] = 0.0;
                }

                for (label jAngle = 0; jAngle < nAngle_; jAngle++)
                {
                    const label rayJ = jAngle + iBand*nAngle_;
                    if (rayJ == rayI)
                    {
                        continue;
                    }
                    const scalar pf = row[jAngle];

                    const scalarField& Ij =
                        ISnapshot_[rayJ].primitiveField();
                    forAll(dsCell, celli)
                    {
                        dsCell[celli] += pf*Ij[celli];
                    }

                    const volScalarField::Boundary& IjBf =
                        ISnapshot_[rayJ].boundaryField();
                    forAll(dsBf, patchi)
                    {
                        scalarField& dsP = dsBf[patchi];
                        const scalarField& IjP = IjBf[patchi];
                        forAll(dsP, fi)
                        {
                            dsP[fi] += pf*IjP[fi];
                        }
                    }
                }
            }

            IRay_[rayI].updateBoundary();
            const scalar maxBandResidual = IRay_[rayI].correct(ds);
            maxResidual = max(maxBandResidual, maxResidual);
        }
    } while (maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}


void Foam::optical::DOM::updateG()
{
    const dimensionedScalar zeroIrradiance("zero", dimPower/dimArea, 0.0);
    G_ = zeroIrradiance;
    forAll(GLambda_, iBand)
    {
        GLambda_[iBand] = zeroIrradiance;
        for (label iAngle = 0; iAngle < nAngle_; iAngle++)
        {
            const label rayI = iAngle + iBand*nAngle_;
            // Convert per-ray radiance [W/m^2/sr] to irradiance [W/m^2]
            // by multiplying by the ray's solid angle.
            GLambda_[iBand] += IRay_[rayI].I()*IRay_[rayI].omega();
        }
        G_ += GLambda_[iBand];
    }
}


const Foam::optical::DOM&
Foam::optical::DOM::lookup(const objectRegistry& db)
{
    return refCast<const DOM>
    (
        db.lookupObject<radiationModel>("opticalRadiationProperties")
    );
}


Foam::label Foam::optical::DOM::nameToRayId(const word& name) const
{
    // assuming name is in the form: CHARS_iBand_iAngle

    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    label ib = readLabel(IStringStream(name.substr(i1+1, i2-i1-1))());
    label ia = readLabel(IStringStream(name.substr(i2+1))());

    return nAngle_*ib + ia;
}


Foam::scalar Foam::optical::DOM::dirToTheta(const vector& dir) const
{
    return Foam::acos(dir.z()/mag(dir));
}


Foam::scalar Foam::optical::DOM::dirToPhi(const vector& dir) const
{
    // atan2 returns phi in (-pi, pi]; map to [0, 2 pi) for compatibility
    // with deltaPhi-based ray binning in dirToRayId. atan2(0, 0) is
    // implementation-defined but conventionally returns 0, which is the
    // correct azimuth for the +z pole (dirToRayId then uses theta to
    // resolve).
    const scalar phi = Foam::atan2(dir.y(), dir.x());
    return phi >= 0.0 ? phi : phi + 2.0*pi;
}


Foam::vector Foam::optical::DOM::anglesToDir
(
    const scalar& theta,
    const scalar& phi
) const
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);
    return vector(sinTheta*cosPhi, sinTheta*sinPhi, cosTheta);
}


Foam::label Foam::optical::DOM::dirToRayId
(
    const vector& dir,
    const label& iBand
) const
{
    scalar tTheta = dirToTheta(dir);
    scalar tPhi = dirToPhi(dir);
    // Clamp to last valid cell so directions at theta=pi or phi=2*pi
    // (e.g. dir = -z, or floating-point round-off near the seam) are
    // assigned to the last angular cell rather than overflowing into
    // the next band.
    label iPhi = min(label(tPhi/deltaPhi_), 2*nPhi_ - 1);
    label iTheta = min(label(tTheta/deltaTheta_), nTheta_ - 1);
    return nAngle_*iBand + iTheta*2*nPhi_ + iPhi;
}


Foam::vector Foam::optical::DOM::intDirOmega
(
    const scalar& theta,
    const scalar& phi,
    const scalar& deltaTheta,
    const scalar& deltaPhi
) const
{
    return vector
        (
            Foam::cos(phi)*Foam::sin(0.5*deltaPhi)
            *(deltaTheta - Foam::cos(2.0*theta)*Foam::sin(deltaTheta)),
            Foam::sin(phi)*Foam::sin(0.5*deltaPhi)
            *(deltaTheta - Foam::cos(2.0*theta)*Foam::sin(deltaTheta)),
            0.5*deltaPhi*Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
        );
}


Foam::vector Foam::optical::DOM::intDirOmega
(
    const scalar& theta,
    const scalar& phi
) const
{
    return intDirOmega(theta, phi, deltaTheta(), deltaPhi());
}


void Foam::optical::DOM::checkDim_()
{
    if (mesh_.nSolutionD() == 2) // 2D (X & Y)
    {
        if (mesh_.solutionD()[vector::Z] != -1)
        {
            FatalErrorInFunction
                << "Currently 2D solution is limited to the x-y plane"
                << exit(FatalError);
        }
        if (nTheta_ != 1)
        {
            FatalErrorInFunction
                << "There must be one theta angle for 2D simulations"
                << exit(FatalError);
        }
    }
    if (mesh_.nSolutionD() == 1)
    {
        FatalErrorInFunction
            << "1D simulations are not supported by DOM; use 2D or 3D"
            << exit(FatalError);
    }
}


void Foam::optical::DOM::setRay_
(
    const label i,
    const label iBand,
    const label iAngle,
    const scalar theta,
    const scalar phi
)
{
    IRay_.set
    (
        i,
        new ray
        (
            *this,
            mesh_,
            iBand,
            iAngle,
            theta,
            phi,
            deltaTheta_,
            deltaPhi_
        )
    );
}

// ************************************************************************* //
