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

Class
    Foam::optical::molecularAbsorptionExtinction

Description
    Generic per-band molecular absorption from a tabulated cross-section
    sigma(lambda_band). The user supplies one pre-band-averaged
    cross-section per radiation band; choice of averaging window
    (line-by-line vs. broadband) is the user's responsibility. Suitable
    for narrow-band UV cases (e.g. KrCl* 222 nm filtered lamps in the
    Herzberg continuum / Hartley band) where the cross-section is
    smooth across the band; less appropriate for broadband UV-C through
    the Schumann-Runge bands without external pre-averaging.

    Two concentration modes:

      concentration idealGas;     uniform N(T,p) * moleFraction
      concentration field;        species field [mol/m^3], spatially varying

    Pure absorption -- scattering is zero. Compose with `rayleigh`
    (for molecular scattering) or another model via `composite`.

    Dictionary:

        extinctionModel molecularAbsorption;
        molecularAbsorptionCoeffs
        {
            nBands           1;
            crossSections    (4.4e-23);   // m^2/molecule, per band

            concentration    idealGas;    // -- or -- field
            T                293.15;      // K           (idealGas only)
            p                101325;      // Pa          (idealGas only)
            moleFraction     0.2095;      // dimensionless  (idealGas only)

            // -- when concentration field --
            field            O3;          // volScalarField in mol/m^3
        }

    Per-band absorption coefficient:

      idealGas:  kappa_b = (p * x / (k_B T)) * sigma_b              [1/m]
      field:     kappa_b = c_field * N_A * sigma_b                  [1/m]

SourceFiles
    molecularAbsorptionExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef opticalRadiationMolecularAbsorptionExtinction_H
#define opticalRadiationMolecularAbsorptionExtinction_H

#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
              Class molecularAbsorptionExtinction Declaration
\*---------------------------------------------------------------------------*/

class molecularAbsorptionExtinction
:
    public extinctionModel
{
public:

    enum class concentrationMode
    {
        idealGas,
        field
    };

private:

    // Private data

        //- Extinction model dictionary
        dictionary coeffsDict_;

        //- Number of bands
        label nBands_;

        //- Concentration mode
        concentrationMode mode_;

        //- Per-band cross-section [m^2/molecule]
        List<scalar> crossSections_;

        //- Gas temperature [K]                      (idealGas mode only)
        scalar T_;

        //- Gas pressure [Pa]                        (idealGas mode only)
        scalar p_;

        //- Mole fraction of absorber in mixture     (idealGas mode only)
        scalar moleFraction_;

        //- Species field name [mol/m^3]             (field mode only)
        word fieldName_;

        //- Owns the species field if not registered (field mode only)
        autoPtr<volScalarField> ownedField_;

public:

    //- Runtime type information
    TypeName("molecularAbsorption");


    // Constructors

        //- Construct from components
        molecularAbsorptionExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~molecularAbsorptionExtinction();


    // Member Functions

        //- Recompute ALambda_ from the species field (no-op in idealGas mode)
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
