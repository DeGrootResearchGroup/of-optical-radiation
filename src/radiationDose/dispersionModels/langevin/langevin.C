/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "langevin.H"
#include "addToRunTimeSelectionTable.H"
#include "gaussianSample.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(langevin, 0);
    addToRunTimeSelectionTable(dispersionModel, langevin, dictionary);

    namespace
    {
        const dispersionModel::addRequiredFields _langevinReqFields
        (
            "langevin",
            &langevin::requiredFields
        );
    }
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::langevin::langevin
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dispersionModel(dict, mesh),
    turbulence_(dict, mesh),
    maxStepFraction_(dict.lookupOrDefault<scalar>("maxStepFraction", 0.02)),
    minLagrangianTime_
    (
        dict.lookupOrDefault<scalar>("minLagrangianTime", 1e-2)
    )
{
    if (maxStepFraction_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "langevin: maxStepFraction must be > 0, got "
            << maxStepFraction_ << exit(FatalIOError);
    }
    if (minLagrangianTime_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "langevin: minLagrangianTime must be >= 0, got "
            << minLagrangianTime_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dose::langevin::integrate
(
    const LangevinState& s,
    scalar dt,
    vector& Xm,
    vector& Xn,
    vector& vEnd
) const
{
    // With the coefficients frozen, v is an Ornstein-Uhlenbeck process
    // relaxing to T grad(sigma) with unit stationary variance. Over dt,
    // h = dt/T, its random part at the end, n, and the random part of its
    // time integral, X, are jointly Gaussian with
    //     var(n)   = 1 - e^-2h
    //     var(X)   = 2 T^2 (h - 2 (1 - e^-h) + (1 - e^-2h)/2)
    //     cov(X,n) = T (1 - e^-h)^2,
    // drawn here as n = sqrt(var n) xi1 and X = (cov/var n) n + r xi2,
    // r^2 = var X - cov^2/var n (series for small h, where the closed form
    // cancels: r^2 = T^2 (h^3/6 - h^5/60)).
    const scalar T = s.TL;
    const scalar h = dt/T;
    const scalar oneMinusE = -expm1(-h);
    const scalar varN = -expm1(-2*h);

    scalar residual;
    if (h < 1e-2)
    {
        residual = sqr(T)*(pow3(h)/6 - pow5(h)/60);
    }
    else
    {
        const scalar varX = 2*sqr(T)*(h - 2*oneMinusE + 0.5*varN);
        const scalar covXn = T*sqr(oneMinusE);
        residual = varX - sqr(covXn)/varN;
    }
    const scalar slope = T*sqr(oneMinusE)/varN;

    const vector n = sqrt(varN)*s.xi1;
    Xn = slope*n + sqrt(max(residual, scalar(0)))*s.xi2;

    // Means: v relaxes from v0 to T grad(sigma)
    const vector vInf = T*s.gradSigma;
    Xm = T*oneMinusE*s.v + (dt - T*oneMinusE)*vInf;
    vEnd = (1 - oneMinusE)*s.v + oneMinusE*vInf + n;

    // The mean displacement that freezing sigma and T_L leaves out, to
    // first order in their gradients with v at its stationary statistics
    // (in units of sigma): sigma varying along the path adds
    // psi1 dt T grad(sigma), and T_L varying along it adds
    // psi2 dt sigma grad(T_L), with
    //     psi1 = 1 - (1 - e^-h)/h,   psi2 = 1 - 2 (1 - e^-h)/h + e^-h.
    // Both tend to one for h >> 1, where the step's mean drift becomes
    // grad(K)/sigma (the random displacement model's), and vanish for
    // h << 1 (series there, where the closed forms cancel).
    // T_L varying along the path also leaves a mean velocity at the step's
    // end, chi sigma grad(T_L) with chi = (1 - e^-h) - h e^-h, which a
    // continuous path carries into the steps that follow; without it the
    // T_L gradient's drift falls short for h ~ 1. sigma varying along the
    // path acts on the displacement only, and leaves none.
    scalar psi1, psi2, chi;
    if (h < 1e-2)
    {
        psi1 = h/2 - sqr(h)/6 + pow3(h)/24;
        psi2 = sqr(h)/6 - pow3(h)/12 + pow4(h)/40;
        chi = sqr(h)/2 - pow3(h)/3 + pow4(h)/8;
    }
    else
    {
        psi1 = 1 - oneMinusE/h;
        psi2 = 1 - 2*oneMinusE/h + (1 - oneMinusE);
        chi = oneMinusE - h*(1 - oneMinusE);
    }
    Xm += dt*(psi1*T*s.gradSigma + psi2*s.sigma*s.gradTL);
    vEnd += chi*s.sigma*s.gradTL;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dose::langevin::correct()
{
    turbulence_.correct();
}


void Foam::dose::langevin::beginStep
(
    dispersionModel::State& state,
    const barycentric& coordinates,
    const tetIndices& tetIs,
    scalar dtMax,
    randomGenerator& rng
) const
{
    LangevinState& s = dynamic_cast<LangevinState&>(state);

    const eddyDiffusivity::sample t = turbulence_.at(coordinates, tetIs);
    s.sigma = t.sigma;
    s.gradSigma = t.gradSigma;
    s.TL = 0.5*t.tau;
    s.gradTL = 0.5*t.gradTau;
    s.K = t.K;

    // A new track takes its velocity from the stationary distribution
    if (!s.hasVelocity)
    {
        s.v = gaussianTriple(rng);
        s.hasVelocity = true;
    }

    s.xi1 = gaussianTriple(rng);
    s.xi2 = gaussianTriple(rng);
}


Foam::vector Foam::dose::langevin::fluctuation
(
    const dispersionModel::State& state,
    scalar dt
) const
{
    const LangevinState& s = dynamic_cast<const LangevinState&>(state);

    vector Xm, Xn, vEnd;
    integrate(s, dt, Xm, Xn, vEnd);
    return s.sigma*(Xm + Xn)/dt;
}


Foam::vector Foam::dose::langevin::cflFluctuation
(
    const dispersionModel::State& state,
    scalar dt
) const
{
    const LangevinState& s = dynamic_cast<const LangevinState&>(state);

    vector Xm, Xn, vEnd;
    integrate(s, dt, Xm, Xn, vEnd);
    return s.sigma*Xm/dt;
}


Foam::scalar Foam::dose::langevin::maxStep
(
    const dispersionModel::State& state,
    scalar maxDisplacement
) const
{
    const LangevinState& s = dynamic_cast<const LangevinState&>(state);

    // Continuous in space (see the class description): no switch, and the
    // T_L floor and the displacement bound are both continuous functions
    // of position.
    const scalar langevinStep =
        maxStepFraction_*max(s.TL, minLagrangianTime_);
    const scalar spreadStep =
        s.K > vSmall ? sqr(maxDisplacement)/(2*s.K) : great;
    return min(langevinStep, spreadStep);
}


void Foam::dose::langevin::endStep
(
    dispersionModel::State& state,
    scalar dt
) const
{
    LangevinState& s = dynamic_cast<LangevinState&>(state);

    vector Xm, Xn, vEnd;
    integrate(s, dt, Xm, Xn, vEnd);
    s.v = vEnd;
}


void Foam::dose::langevin::reflect
(
    dispersionModel::State& state,
    const vector& n
) const
{
    LangevinState& s = dynamic_cast<LangevinState&>(state);
    s.v -= 2.0*(s.v & n)*n;
}


// ************************************************************************* //
