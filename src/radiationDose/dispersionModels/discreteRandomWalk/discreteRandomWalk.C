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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::discreteRandomWalk::discreteRandomWalk
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dispersionModel(dict, mesh),
    turbulence_(dict, mesh),
    wellMixed_(dict.lookupOrDefault<Switch>("wellMixed", true))
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dose::discreteRandomWalk::beginStepUncorrected
(
    DRWState& s,
    label celli,
    scalar dtMax,
    randomGenerator& rng
) const
{
    scalar kVal, tauE;
    turbulence_.cell(celli, kVal, tauE);

    if (s.remaining <= 0 || kVal <= small)
    {
        // Resample: 3 independent N(0, sigma) components from a shared
        // Box-Muller triple. sigma = sqrt(2k/3) is the per-component
        // isotropic-turbulence variance (Gosman-Ioannides 1981).
        const scalar sigma = sqrt(2.0/3.0*kVal);
        s.uPrime = sigma*gaussianTriple(rng);
        s.remaining = tauE;
    }

    s.remaining -= dtMax;
}


bool Foam::dose::discreteRandomWalk::newEddies
(
    const DRWState& s,
    scalar dt,
    scalar& m,
    scalar& frac
) const
{
    if (s.sigma <= 0 || s.tau <= vSmall)
    {
        return false;
    }

    const scalar tLeft = dt - max(s.remaining, scalar(0));
    m = floor(tLeft/s.tau);
    frac = max(tLeft - m*s.tau, scalar(0));
    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dose::discreteRandomWalk::correct()
{
    if (wellMixed_)
    {
        turbulence_.correct();
    }
}


void Foam::dose::discreteRandomWalk::beginStep
(
    dispersionModel::State& state,
    const barycentric& coordinates,
    const tetIndices& tetIs,
    scalar dtMax,
    randomGenerator& rng
) const
{
    // The track owns its own DRWState. No locking needed because each
    // track is integrated by exactly one thread at a time.
    DRWState& s = dynamic_cast<DRWState&>(state);

    if (!wellMixed_)
    {
        beginStepUncorrected(s, tetIs.cell(), dtMax, rng);
        return;
    }

    const eddyDiffusivity::sample t = turbulence_.at(coordinates, tetIs);
    s.sigma = t.sigma;
    s.tau = t.tau;
    s.drift = t.gradK;

    // A step of up to dtMax may outlast the current eddy: draw the sum
    // of the whole eddies that fit and the eddy that carries over, unless
    // an earlier step drew them and did not use them. Isotropic N(0, 1)
    // triples, scaled in fluctuation().
    if (s.remaining < dtMax && !s.drawn)
    {
        s.xiSum = gaussianTriple(rng);
        s.xiNew = gaussianTriple(rng);
        s.drawn = true;
    }
}


Foam::vector Foam::dose::discreteRandomWalk::fluctuation
(
    const dispersionModel::State& state,
    scalar dt
) const
{
    const DRWState& s = dynamic_cast<const DRWState&>(state);

    if (!wellMixed_)
    {
        return s.uPrime;
    }

    if (s.remaining >= dt)
    {
        return s.uPrime + s.drift;
    }

    // The step outlasts the current eddy. Its displacement is the rest
    // of that eddy, then m whole eddies of duration tau (independent
    // N(0, sigma^2) velocities, so their summed displacement is one
    // Gaussian of standard deviation sigma tau sqrt(m)), then the first
    // `frac` of the eddy that carries into the next step.
    vector disp = max(s.remaining, scalar(0))*s.uPrime;

    scalar m = 0, frac = 0;
    if (newEddies(s, dt, m, frac))
    {
        disp += s.sigma*(s.tau*sqrt(m)*s.xiSum + frac*s.xiNew);
    }

    return disp/dt + s.drift;
}


void Foam::dose::discreteRandomWalk::endStep
(
    dispersionModel::State& state,
    scalar dt
) const
{
    if (!wellMixed_)
    {
        return;
    }

    DRWState& s = dynamic_cast<DRWState&>(state);

    if (s.remaining >= dt)
    {
        s.remaining -= dt;
        return;
    }

    // The step used the draws: the eddy that carries over begins
    scalar m = 0, frac = 0;
    if (newEddies(s, dt, m, frac))
    {
        s.uPrime = s.sigma*s.xiNew;
        s.remaining = s.tau - frac;
    }
    else
    {
        s.uPrime = vector::zero;
        s.remaining = 0;
    }
    s.drawn = false;
}


void Foam::dose::discreteRandomWalk::reflect
(
    dispersionModel::State& state,
    const vector& n
) const
{
    if (!wellMixed_)
    {
        return;
    }

    // The eddy the particle carries, and the one that would carry over
    // from this step, are mirrored with the particle's own velocity.
    DRWState& s = dynamic_cast<DRWState&>(state);
    s.uPrime -= 2.0*(s.uPrime & n)*n;
    s.xiNew -= 2.0*(s.xiNew & n)*n;
}


// ************************************************************************* //
