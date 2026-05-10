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

#include "molecularAbsorptionExtinction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(molecularAbsorptionExtinction, 0);

        addToRunTimeSelectionTable
        (
            extinctionModel,
            molecularAbsorptionExtinction,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::molecularAbsorptionExtinction::molecularAbsorptionExtinction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    extinctionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    mode_(concentrationMode::idealGas),
    crossSections_(),
    T_(0.0),
    p_(0.0),
    moleFraction_(0.0),
    fieldName_(),
    ownedField_()
{
    // init() sets the inherited nBands_ and resizes ALambda_ / SLambda_.
    init(readLabel(coeffsDict_.lookup("nBands")));

    coeffsDict_.lookup("crossSections") >> crossSections_;

    if (crossSections_.size() != nBands_)
    {
        FatalErrorIn
        (
            "molecularAbsorptionExtinction::"
            "molecularAbsorptionExtinction(const dictionary&, const fvMesh&)"
        )   << "crossSections list has " << crossSections_.size()
            << " entries but nBands = " << nBands_
            << exit(FatalError);
    }

    const word concKw(coeffsDict_.lookup("concentration"));

    if (concKw == "idealGas")
    {
        mode_ = concentrationMode::idealGas;
        T_            = readScalar(coeffsDict_.lookup("T"));
        p_            = readScalar(coeffsDict_.lookup("p"));
        moleFraction_ = readScalar(coeffsDict_.lookup("moleFraction"));

        const scalar kB = 1.380649e-23;
        const scalar N  = p_*moleFraction_/(kB*T_);

        forAll(crossSections_, b)
        {
            const scalar kappaB = N*crossSections_[b];
            ALambda_[b] = dimensionedScalar("A", dimless/dimLength, kappaB);

            Info<< "    molecularAbsorption(idealGas) band " << b
                << ": sigma = " << crossSections_[b]
                << " m^2/mol., N = " << N << " 1/m^3, "
                << "kappa = " << kappaB << " 1/m" << endl;
        }
    }
    else if (concKw == "field")
    {
        mode_ = concentrationMode::field;
        fieldName_ = word(coeffsDict_.lookup("field"));

        // Mirror linearSpeciesExtinction: load from disk only when an
        // upstream solver hasn't already registered the field.
        if (!mesh.foundObject<volScalarField>(fieldName_))
        {
            ownedField_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        fieldName_,
                        mesh.time().name(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                )
            );
        }

        forAll(crossSections_, b)
        {
            Info<< "    molecularAbsorption(field=" << fieldName_
                << ") band " << b
                << ": sigma = " << crossSections_[b] << " m^2/mol." << endl;
        }

        correct();
    }
    else
    {
        FatalErrorIn
        (
            "molecularAbsorptionExtinction::"
            "molecularAbsorptionExtinction(const dictionary&, const fvMesh&)"
        )   << "concentration must be 'idealGas' or 'field', got '"
            << concKw << "'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::molecularAbsorptionExtinction::~molecularAbsorptionExtinction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optical::molecularAbsorptionExtinction::correct()
{
    if (mode_ == concentrationMode::idealGas)
    {
        // ALambda_ was set in the constructor and is uniform/constant.
        return;
    }

    const scalar Na = 6.02214076e23;

    const volScalarField& cField =
        mesh_.lookupObject<volScalarField>(fieldName_);

    forAll(ALambda_, b)
    {
        const dimensionedScalar coeff
        (
            "sigmaNa",
            dimensionSet(0, 2, 0, 0, -1, 0, 0),  // m^2/mol
            Na*crossSections_[b]
        );

        ALambda_[b] = coeff*cField;
    }
}


// ************************************************************************* //
