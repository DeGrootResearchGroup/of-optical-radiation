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

#include "linearSpeciesExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(linearSpeciesExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            linearSpeciesExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::linearSpeciesExtinction::linearSpeciesExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    absorption_(coeffsDict_.lookup("absorption")),
    nAbsorbing_(readLabel(coeffsDict_.lookup("nAbsorbing"))),
    aSpecies_(nAbsorbing_),
    aCoeffs_(nAbsorbing_),
    scattering_(coeffsDict_.lookup("scattering")),
    nScattering_(readLabel(coeffsDict_.lookup("nScattering"))),
    sSpecies_(nScattering_),
    sCoeffs_(nScattering_)
{
    // init() sets the inherited nBands_ and resizes ALambda_ / SLambda_.
    init(readLabel(coeffsDict_.lookup("nBands")));

    forAll(aCoeffs_, i)
    {
        aCoeffs_[i].setSize(nBands_, 0.0);
    }
    forAll(sCoeffs_, i)
    {
        sCoeffs_[i].setSize(nBands_, 0.0);
    }

    // Read the coefficients
    if (absorption_) coeffsDict_.lookup("absorptionCoeff") >> aCoeffs_;
    if (scattering_) coeffsDict_.lookup("scatteringCoeff") >> sCoeffs_;

    // Read the species names
    if (absorption_) coeffsDict_.lookup("absorbingVars") >> aSpecies_;
    if (scattering_) coeffsDict_.lookup("scatteringVars") >> sSpecies_;

    // Load any species fields that aren't already registered with the
    // mesh (e.g. when running standalone without a coupled species
    // solver). Fields constructed here auto-register with the mesh
    // registry, so correct()'s lookupObject calls find them. We hold
    // ownership in ownedSpeciesFields_ so they live as long as we do.
    DynamicList<word> speciesToLoad;
    if (absorption_)
    {
        forAll(aSpecies_, i)
        {
            if (!mesh.foundObject<volScalarField>(aSpecies_[i]))
            {
                speciesToLoad.append(aSpecies_[i]);
            }
        }
    }
    if (scattering_)
    {
        forAll(sSpecies_, i)
        {
            if (!mesh.foundObject<volScalarField>(sSpecies_[i]))
            {
                speciesToLoad.append(sSpecies_[i]);
            }
        }
    }

    ownedSpeciesFields_.setSize(speciesToLoad.size());
    forAll(speciesToLoad, i)
    {
        ownedSpeciesFields_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    speciesToLoad[i],
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
    }

    // Validate species-field dimensions at construction so dictionary
    // mistakes (e.g. supplying a mole-fraction field instead of a mass
    // concentration) are caught upfront rather than at the first
    // correct() call, deep inside the time loop. Expected dimensions:
    // mass concentration [kg/m^3] so that
    //   a * C = [m^2/kg] * [kg/m^3] = [1/m]
    // matches the extinction-coefficient dimensions assigned to
    // ALambda_ / SLambda_. (See the dimensionedScalar `a`/`s`
    // constructions in correct() below.)
    const dimensionSet expectedDims(dimMass/pow3(dimLength));
    auto checkDim = [&](const word& name)
    {
        const volScalarField& f = mesh.lookupObject<volScalarField>(name);
        if (f.dimensions() != expectedDims)
        {
            FatalErrorInFunction
                << "Species field '" << name << "' has dimensions "
                << f.dimensions() << "; linearSpeciesExtinction expects "
                << "mass-concentration dimensions " << expectedDims
                << " (kg/m^3). Coefficient * concentration must yield"
                << " an extinction coefficient with 1/m dimensions."
                << exit(FatalError);
        }
    };
    if (absorption_)
    {
        forAll(aSpecies_, i) checkDim(aSpecies_[i]);
    }
    if (scattering_)
    {
        forAll(sSpecies_, i) checkDim(sSpecies_[i]);
    }

    // Correct the extinction coefficient fields
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::optical::linearSpeciesExtinction::~linearSpeciesExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::optical::linearSpeciesExtinction::correct()
{
    // Absorption / scattering loops are guarded on the master toggle.
    // If absorption is off the user is allowed to leave absorbingVars
    // / absorptionCoeff unset; without this guard the loop would still
    // run nAbsorbing_ times (whatever the user wrote in the dict) and
    // try to look up empty species names. Symmetric for scattering.
    forAll(ALambda_, iBand)
    {
        ALambda_[iBand] = dimensionedScalar("A", dimless/dimLength, 0.0);
        if (!absorption_)
        {
            continue;
        }
        for (label i = 0; i < nAbsorbing_; i++)
        {
            const dimensionedScalar a
            (
                "a", dimLength*dimLength/dimMass, aCoeffs_[i][iBand]
            );
            ALambda_[iBand] +=
                a*mesh().lookupObject<volScalarField>(aSpecies_[i]);
        }
    }

    forAll(SLambda_, iBand)
    {
        SLambda_[iBand] = dimensionedScalar("S", dimless/dimLength, 0.0);
        if (!scattering_)
        {
            continue;
        }
        for (label i = 0; i < nScattering_; i++)
        {
            const dimensionedScalar s
            (
                "s", dimLength*dimLength/dimMass, sCoeffs_[i][iBand]
            );
            SLambda_[iBand] +=
                s*mesh().lookupObject<volScalarField>(sSpecies_[i]);
        }
    }
}


// ************************************************************************* //


