/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Class
    Foam::optical::mieExtinction

Description
    Mie scattering extinction model for monodisperse spherical particles
    suspended in a host medium. Per-band absorption and scattering
    coefficients are

        kappa(x)   = pi r^2 Q_abs(lambda_b) N(x)
        sigma_s(x) = pi r^2 Q_sca(lambda_b) N(x)

    where r is the particle radius, Q_abs and Q_sca are Mie
    absorption / scattering efficiencies (computed once at construction
    via mieKernel for size parameter
    x_b = 2 pi r n_medium / lambda_b and complex relative refractive
    index m_rel = m_particle / n_medium), and N(x) is the local particle
    number density read from a registered volScalarField (units 1/m^3).

    The number-density field is loaded from disk if no other component
    has registered it (mirrors linearSpeciesExtinction).

    Limitations (see CLAUDE.md for context):
      - Monodisperse: a single radius is used for the whole field. A
        polydisperse extension would integrate the efficiencies over
        a size distribution at construction.
      - Mass / volume concentration must be converted to number density
        by the user (e.g. n = c / (rho_p * 4/3 pi r^3) for a mass
        concentration c [kg/m^3]).

SourceFiles
    mieExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef optical_mieExtinction_H
#define optical_mieExtinction_H

#include "extinctionModel.H"

#include <complex>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                       Class mieExtinction Declaration
\*---------------------------------------------------------------------------*/

class mieExtinction
:
    public extinctionModel
{
    // Private data

        //- Extinction model dictionary
        dictionary coeffsDict_;

        //- Per-band central wavelengths in vacuum [nm]
        List<scalar> wavelengths_;

        //- Particle radius [m]
        scalar radius_;

        //- Complex particle refractive index (real, imag)
        std::complex<scalar> mParticle_;

        //- Real refractive index of the host medium
        scalar mMedium_;

        //- Name of the volScalarField holding particle number density [1/m^3]
        word numberDensityFieldName_;

        //- Number-density field loaded by us (held for lifetime if not
        //  registered by an external solver)
        autoPtr<volScalarField> ownedNumberDensityField_;

        //- Scattering efficiency Q_sca per band (cached for write())
        List<scalar> QscaPerBand_;

        //- Absorption efficiency Q_abs per band (cached for write())
        List<scalar> QabsPerBand_;

        //- Pre-multiplied geometric factor pi r^2 Q_sca per band [m^2]
        List<scalar> sigmaScaPrefactorPerBand_;

        //- Pre-multiplied geometric factor pi r^2 Q_abs per band [m^2]
        List<scalar> sigmaAbsPrefactorPerBand_;

public:

    //- Runtime type information
    TypeName("mie");


    // Constructors

        //- Construct from components
        mieExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~mieExtinction();


    // Member Functions

        //- Update kappa and sigma_s fields from the current number density
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
