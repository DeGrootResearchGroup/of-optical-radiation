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

\*---------------------------------------------------------------------------*/

#include "rayleighExtinction.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(rayleighExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            rayleighExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * * //

namespace
{
    // Peck & Reeder 1972 dispersion fit for dry air at 288.15 K, 101325 Pa
    // (n - 1) form, with sigma the wavenumber [1/um]. The fit's quoted range
    // is 230-1690 nm; extrapolation to 222 nm shifts (n - 1) by less than
    // 2 %, well within the (n^2 - 1)^2 sensitivity.
    inline Foam::scalar airRefractiveIndex(const Foam::scalar lambdaNm)
    {
        const Foam::scalar sigma2 = 1.0e6/(lambdaNm*lambdaNm); // 1/um^2
        return
            1.0
          + 8060.51e-8
          + 2480990.0e-8/(132.274 - sigma2)
          + 17455.7e-8/(39.32957 - sigma2);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::rayleighExtinction::rayleighExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    wavelengths_(),
    T_(coeffsDict_.lookupOrDefault<scalar>("T", 288.15)),
    p_(coeffsDict_.lookupOrDefault<scalar>("p", 101325.0)),
    kingFactor_(coeffsDict_.lookupOrDefault<scalar>("kingFactor", 1.05)),
    sigmaSPerBand_()
{
    // init() sets the inherited nBands_ and resizes ALambda_ / SLambda_.
    init(readLabel(coeffsDict_.lookup("nBands")));
    sigmaSPerBand_.setSize(nBands_, 0.0);

    coeffsDict_.lookup("wavelengths") >> wavelengths_;

    if (wavelengths_.size() != nBands_)
    {
        FatalErrorIn
        (
            "rayleighExtinction::rayleighExtinction"
            "(const dictionary&, const fvMesh&)"
        )   << "wavelengths list has " << wavelengths_.size()
            << " entries but nBands = " << nBands_
            << exit(FatalError);
    }

    const scalar kB = 1.380649e-23;
    const scalar pi = constant::mathematical::pi;

    // Bodhaine et al. 1999 reference state: the standard conditions at
    // which the Peck-Reeder air refractive-index fit (used in
    // airRefractiveIndex above) is defined. The Rayleigh per-molecule
    // cross-section is intrinsic to the molecule and was inferred from
    // n_s measured at N_s, so the (n^2 - 1) and N^2 factors in
    // sigma(lambda) must use that reference density -- NOT the local
    // density implied by the user-supplied (T, p). Local thermodynamics
    // enters only through the bulk scattering coefficient
    //   sigma_s(T, p) = N(T, p) * sigma_per_molecule
    // below. Using local N in the cross-section formula gives results
    // that are correct at the Peck-Reeder reference state but drift
    // (~few %) away from it.
    const scalar T_s = 288.15;       // Peck & Reeder reference T [K]
    const scalar p_s = 101325.0;     // Peck & Reeder reference p [Pa]
    const scalar N_s = p_s/(kB*T_s); // standard-state number density [1/m^3]
    const scalar N   = p_/(kB*T_);   // local number density [1/m^3]

    forAll(wavelengths_, b)
    {
        const scalar lamM = wavelengths_[b]*1e-9;
        const scalar ns   = airRefractiveIndex(wavelengths_[b]);
        const scalar lr   = (sqr(ns) - 1.0)/(sqr(ns) + 2.0);
        // Per-molecule cross-section (Bodhaine 1999), uses N_s.
        const scalar sigmaPerMolecule =
            24.0*pow3(pi)/(pow4(lamM)*sqr(N_s))*sqr(lr)*kingFactor_;

        // Bulk volume scattering coefficient at local (T, p).
        sigmaSPerBand_[b] = N*sigmaPerMolecule;

        ALambda_[b] = dimensionedScalar("A", dimless/dimLength, 0.0);
        SLambda_[b] =
            dimensionedScalar("S", dimless/dimLength, sigmaSPerBand_[b]);

        Info<< "    rayleigh band " << b
            << ": lambda = " << wavelengths_[b] << " nm, "
            << "sigma_s = " << sigmaSPerBand_[b] << " 1/m" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::rayleighExtinction::~rayleighExtinction()
{}


// ************************************************************************* //
