/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "discreteRandomWalk.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(discreteRandomWalk, 0);
    addToRunTimeSelectionTable(dispersionModel, discreteRandomWalk, dictionary);
}
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
    epsilonName_(dict.lookupOrDefault<word>("epsilon", "epsilon")),
    Cl_(dict.lookupOrDefault<scalar>("Cl", 0.15))
{}


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
    const volScalarField& epsField =
        mesh_.lookupObject<volScalarField>(epsilonName_);

    const scalar kVal = max(kField[celli], scalar(0));
    const scalar epsVal = max(epsField[celli], small);

    // The track owns its own DRWState. No locking needed because each
    // track is integrated by exactly one thread at a time.
    DRWState& s = dynamic_cast<DRWState&>(state);

    if (s.remaining <= 0 || kVal <= small)
    {
        // Resample: 3 independent N(0, sigma) components via Box-Muller.
        const scalar sigma = sqrt(2.0/3.0*kVal);

        const scalar twoPi = constant::mathematical::twoPi;
        const scalar r1 = max(rng.scalar01(), small);
        const scalar r2 = rng.scalar01();
        const scalar r3 = max(rng.scalar01(), small);
        const scalar r4 = rng.scalar01();

        const scalar mag1 = sqrt(-2.0*log(r1));
        const scalar z1 = mag1*cos(twoPi*r2);
        const scalar z2 = mag1*sin(twoPi*r2);
        const scalar mag3 = sqrt(-2.0*log(r3));
        const scalar z3 = mag3*cos(twoPi*r4);

        s.uPrime = sigma*vector(z1, z2, z3);
        s.remaining = Cl_*kVal/epsVal;
    }

    s.remaining -= dt;
    return s.uPrime;
}


// ************************************************************************* //
