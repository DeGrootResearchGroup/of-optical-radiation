/*---------------------------------------------------------------------------*\
  =========   	              |
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

#include "photoBioDOM.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "phaseFunctionModel.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(photoBioDOM, 0);
        addToRunTimeSelectionTable
        (
            photoBioModel,
            photoBioDOM,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::photoBioDOM::photoBioDOM(const volScalarField& intensity)
:
    photoBioModel(typeName, intensity),
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
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    diffusionScatter_
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
        dimensionedScalar("diffusionScatter",dimMass/pow3(dimTime), 0.0)
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
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50))
{
    Info<< "Creating photoBioDOM model with " << nBand_ << " bands" << endl;

    // Check that dimension of mesh is compatible with settings
    checkDim_();

    // Set the number of angles (=2*nPhi*nTheta for 2D/3D, =2 for 1D)
    nAngle_ = mesh_.nSolutionD() > 1 ? 2*nPhi_*nTheta_ : 2;

    // Set the number of rays, and allocate pointer list
    nRay_ = nAngle_*nBand_;
    IRay_.setSize(nRay_);

    // Set deltaTheta (=pi/nTheta for 3D; =pi for 1D/2D)
    deltaTheta_ = mesh_.nSolutionD() == 3 ? pi/nTheta_ : pi;

    // Set deltaPhi (=pi/nPhi for 2D/3D, =pi for 1D)
    deltaPhi_ = mesh_.nSolutionD() > 1 ? pi/nPhi_ : pi;

    // Set up all of the rays
    label i = 0;
    for (label iBand = 0; iBand < nBand_; iBand++)
    {
        scalar sumOmega = 0.0;
        for (label iTheta = 0; iTheta < nTheta_; iTheta++)
        {
            for (label iPhi = 0; iPhi < nAngle_/nTheta_; iPhi++)
            {
                label iAngle = iPhi + 2*iTheta*nPhi_;
                scalar theta = (iTheta + 0.5)*deltaTheta_;
                scalar phi = (iPhi + 0.5)*deltaPhi_;
                setRay_(i, iBand, iAngle, phi, theta);
                sumOmega += IRay_[i].omega();
                i++;
            }
        }
        if (iBand == 0)
        {
            Info << "Sum of solid angles: " << sumOmega << endl;
        }
    }

    Info<< "photoBioDOM : Allocated " << IRay_.size() << " rays" << endl;

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

    if (phaseFunctionModel_->inScatter())
    {
        pf0_.setSize(nBand_ * nAngle_);
        for (label iBand = 0 ; iBand < nBand_; iBand++)
	{
            for (label iAngle = 0 ; iAngle < nAngle_; iAngle++)
            {
                pf0_[iAngle+iBand*nAngle_] = phaseFunctionModel_->correct(iAngle,iAngle,iBand);
            }
        }
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::photoBio::photoBioDOM::~photoBioDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::photoBio::photoBioDOM::read()
{
    if (photoBioModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);
        return true;
    }
    else
    {
        return false;
    }
}

void Foam::photoBio::photoBioDOM::calculate()
{
    scalar maxResidual = 0.0;
    label rayJ;
    label iBand = 0;
    label iAngle = 0;
    label radIter = 0;
    scalar maxBandResidual = 0.0;

    // Correct the extinction model
    extinction_->correct();

    do
    {
        radIter++;
        maxResidual = 0.0;
        forAll(IRay_, rayI)
	{
            iBand = IRay_[rayI].iBand();
            iAngle = IRay_[rayI].iAngle();

            Info << endl;
            Info << "photoBio solver: " << "    iter: " << radIter
                                        << "    iBand: "<< iBand
                                        << "    iAngle  : "<< iAngle;
            Info << endl;

            if (phaseFunctionModel_->inScatter())
            {
                diffusionScatter_ = dimensionedScalar("diffusionScatter",dimMass/pow3(dimTime), 0.0) ;

                for (label jAngle = 0; jAngle < nAngle_; jAngle++)
                {
                    rayJ = jAngle + iBand*nAngle_;
                    if(rayI != rayJ )
                    {
                        scalar rayCos = IRay_[rayI].d() & IRay_[rayJ].d();
                        if(rayCos > 0 )
                        {
                            diffusionScatter_ = diffusionScatter_ + IRay_[rayJ].I()*phaseFunctionModel_->correct(rayI,rayJ, iBand)*IRay_[rayJ].omega();
                        }
                    }
                }
            }

            IRay_[rayI].updateBoundary();
            maxBandResidual = IRay_[rayI].correct(); // solve the RTE equation
            maxResidual = max(maxBandResidual, maxResidual);
        }
    } while(maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}


void Foam::photoBio::photoBioDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    label rayI;
    forAll(GLambda_, iBand)
    {
        GLambda_[iBand] = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
        for (label iAngle = 0; iAngle < nAngle_; iAngle++)
        {
            rayI = iAngle + iBand*nAngle_;
            GLambda_[iBand] += IRay_[rayI].I()*IRay_[rayI].omega();  // convert the radiance to irradiance
        }
        G_ += GLambda_[iBand];
    }
}


Foam::label Foam::photoBio::photoBioDOM::nameToRayId(const word& name) const
{
    // assuming name is in the form: CHARS_iBand_iAngle

    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    label ib = readLabel(IStringStream(name.substr(i1+1, i2-1))());
    label ia = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());

    return nAngle_*ib + ia;
}


Foam::scalar Foam::photoBio::photoBioDOM::dirToTheta(const vector& dir) const
{
    return Foam::acos(dir.z()/mag(dir));
}


Foam::scalar Foam::photoBio::photoBioDOM::dirToPhi(const vector& dir) const
{
    scalar phi = 0.0;
    if (mag(dir.x()) > SMALL)
    {
        phi = Foam::atan(dir.y()/dir.x());
        if (dir.x() < 0.0) phi += pi;
        if (dir.x() > 0.0 && dir.y() < 0.0) phi += 2.0*pi;
    }
    else
    {
        if (dir.y() > 0.0) phi = pi/2.0;
        if (dir.y() < 0.0) phi = 3.0*pi/2.0;
    }
    return phi;
}


Foam::vector Foam::photoBio::photoBioDOM::anglesToDir
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


Foam::label Foam::photoBio::photoBioDOM::dirToRayId
(
    const vector& dir,
    const label& iBand
) const
{
    scalar tTheta = dirToTheta(dir);
    scalar tPhi = dirToPhi(dir);
    label iPhi = label(tPhi/deltaPhi_);
    label iTheta = label(tTheta/deltaTheta_);
    return nAngle_*iBand + iTheta*2*nPhi_ + iPhi;
}


Foam::vector Foam::photoBio::photoBioDOM::intDirOmega
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


Foam::vector Foam::photoBio::photoBioDOM::intDirOmega
(
    const scalar& theta,
    const scalar& phi
) const
{
    return intDirOmega(theta, phi, deltaTheta(), deltaPhi());
}


void Foam::photoBio::photoBioDOM::checkDim_()
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
    if (mesh_.nSolutionD() == 1) // 1D (X)
    {
        if (mesh_.solutionD()[vector::X] != 1)
        {
            FatalErrorInFunction
                << "Currently 1D solution is limited to the x-direction"
                << exit(FatalError);
        }
        if (nTheta_ != 1)
        {
            FatalErrorInFunction
                << "There must be one theta angle for 1D simulations"
                << exit(FatalError);
        }
        if (nPhi_ != 2)
        {
            FatalErrorInFunction
                << "There must be two phi angles for 1D simulations"
                << exit(FatalError);
        }
    }
}


void Foam::photoBio::photoBioDOM::setRay_
(
    const label i,
    const label iBand,
    const label iAngle,
    const scalar phi,
    const scalar theta
)
{
    IRay_.set
    (
        i,
        new photoBioIntensityRay
        (
            *this,
            mesh_,
            iBand,
            iAngle,
            phi,
            theta,
            deltaPhi_,
            deltaTheta_
        )
    );
}

// ************************************************************************* //
