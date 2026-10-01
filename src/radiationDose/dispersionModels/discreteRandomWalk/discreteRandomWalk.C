/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "discreteRandomWalk.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "volFields.H"
#include "gaussianSample.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(discreteRandomWalk, 0);
    addToRunTimeSelectionTable(dispersionModel, discreteRandomWalk, dictionary);

    // Register the static requiredFields helper so the base class's
    // dispatcher (dispersionModel::requiredFields(dict)) can find it
    // without instantiating a throwaway DRW model.
    namespace
    {
        const dispersionModel::addRequiredFields _drwReqFields
        (
            "discreteRandomWalk",
            &discreteRandomWalk::requiredFields
        );
    }
}
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::dose::discreteRandomWalk::dissipationName
(
    const dictionary& dict
)
{
    if (dict.found("omega") && dict.found("epsilon"))
    {
        FatalIOErrorInFunction(dict)
            << "discreteRandomWalk: name either epsilon or omega, not both"
            << exit(FatalIOError);
    }
    return dict.found("omega")
        ? dict.lookup<word>("omega")
        : dict.lookupOrDefault<word>("epsilon", "epsilon");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::discreteRandomWalk::discreteRandomWalk
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dispersionModel(dict, mesh),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    omega_(dict.found("omega")),
    dissipationName_(dissipationName(dict)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Cl_(dict.lookupOrDefault<scalar>("Cl", 0.15)),
    tauEMax_(dict.lookupOrDefault<scalar>("tauEMax", 100.0))
{
    if (Cmu_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "discreteRandomWalk: Cmu must be > 0, got " << Cmu_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::dose::discreteRandomWalk::fluctuation
(
    dispersionModel::State& state,
    const vector& x,
    label celli,
    scalar dt,
    randomGenerator& rng
) const
{
    if (celli < 0)
    {
        return vector::zero;
    }

    const volScalarField& kField =
        mesh_.lookupObject<volScalarField>(kName_);
    const volScalarField& dissField =
        mesh_.lookupObject<volScalarField>(dissipationName_);

    const scalar kVal = max(kField[celli], scalar(0));
    const scalar dissVal = max(dissField[celli], small);

    // The track owns its own DRWState. No locking needed because each
    // track is integrated by exactly one thread at a time.
    DRWState& s = dynamic_cast<DRWState&>(state);

    if (s.remaining <= 0 || kVal <= small)
    {
        // Resample: 3 independent N(0, sigma) components from a shared
        // Box-Muller triple. sigma = sqrt(2k/3) is the per-component
        // isotropic-turbulence variance (Gosman-Ioannides 1981).
        const scalar sigma = sqrt(2.0/3.0*kVal);
        s.uPrime = sigma*gaussianTriple(rng);
        // Cl k / epsilon, which with epsilon = Cmu k omega is Cl / (Cmu omega).
        const scalar tauE = omega_ ? Cl_/(Cmu_*dissVal) : Cl_*kVal/dissVal;
        s.remaining = min(tauE, tauEMax_);
    }

    s.remaining -= dt;
    return s.uPrime;
}


// ************************************************************************* //
