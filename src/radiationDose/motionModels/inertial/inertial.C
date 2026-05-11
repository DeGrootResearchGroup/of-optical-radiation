/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "inertial.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(inertial, 0);
    addToRunTimeSelectionTable(motionModel, inertial, dictionary);

    namespace
    {
        const motionModel::addRequiredFields _inertialReqFields
        (
            "inertial",
            &inertial::requiredFields
        );
    }
}
}


// Boltzmann constant [J/K]. Hand-rolled to keep the dependency surface
// minimal — Foam::constant::physicoChemical exposes it but pulling that
// in drags more than this single number is worth.
namespace
{
    static constexpr Foam::scalar kB = 1.380649e-23;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::inertial::inertial
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    motionModel(dict, mesh),
    drag_(dragModel::New(dict.subDict("drag"))),
    rhoP_(readScalar(dict.lookup("rhoP"))),
    dP_(readScalar(dict.lookup("dP"))),
    rhoF_(dict.lookupOrDefault<scalar>("rhoF", 1000.0)),
    muF_(dict.lookupOrDefault<scalar>("muF", 1.0e-3)),
    nuF_(muF_/rhoF_),
    mP_((constant::mathematical::pi/6.0)*rhoP_*pow3(dP_)),
    gravityActive_(false),
    gravity_(Zero),
    brownianActive_(false),
    T_(293.15)
{
    if (dict.found("gravity"))
    {
        const dictionary& gDict = dict.subDict("gravity");
        gravityActive_ = gDict.lookupOrDefault<Switch>("active", true);
        if (gravityActive_)
        {
            gravity_ = vector(gDict.lookup("value"));
        }
    }

    if (dict.found("brownian"))
    {
        const dictionary& bDict = dict.subDict("brownian");
        brownianActive_ = bDict.lookupOrDefault<Switch>("active", true);
        if (brownianActive_)
        {
            T_ = bDict.lookupOrDefault<scalar>("T", 293.15);
        }
    }

    // Diagnostics: report Stokes settling and a Brownian-only diffusion
    // estimate at construction so dictionary mistakes show up in the
    // log rather than as quiet wrong-answers downstream.
    const scalar tauStokes = rhoP_*dP_*dP_/(18.0*muF_);
    Info<< "  inertial motion: rhoP=" << rhoP_ << " kg/m^3"
        << ", dP=" << dP_ << " m"
        << ", mP=" << mP_ << " kg"
        << ", tau_p (Stokes)=" << tauStokes << " s" << endl;

    if (gravityActive_)
    {
        const vector Vsettling = tauStokes*(1.0 - rhoF_/rhoP_)*gravity_;
        Info<< "  Stokes settling velocity = " << Vsettling
            << " m/s" << endl;
    }

    if (brownianActive_)
    {
        const scalar D = kB*T_*tauStokes/mP_;
        Info<< "  Brownian D (Einstein-Stokes) = " << D
            << " m^2/s at T=" << T_ << " K" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dose::motionModel::stepResult
Foam::dose::inertial::advance
(
    State& /*state*/,
    const vector& V,
    const vector& Umean,
    const vector& uPrime,
    scalar dt,
    const vector& xi
) const
{
    const vector Useen = Umean + uPrime;
    const vector ag =
        gravityActive_
      ? (1.0 - rhoF_/rhoP_)*gravity_
      : vector::zero;

    const scalar tauP = drag_->tauP(V, Useen, rhoP_, dP_, nuF_, muF_);
    const vector Veq  = Useen + ag*tauP;

    const scalar omega = dt/tauP;
    const scalar expw  = exp(-omega);

    // (1 - exp(-omega))/omega; well-conditioned for small omega via
    // expm1, recovers 1 in the omega→0 limit.
    const scalar phi =
        (omega > small)
      ? -expm1(-omega)/omega
      : 1.0 - 0.5*omega;

    // Drag-only end-of-step velocity (OU exact, deterministic part)
    const vector VdragNew = Veq + (V - Veq)*expw;

    // Mean velocity over the step (for displacement)
    const vector Vdisp = Veq + (V - Veq)*phi;

    // Brownian kick on V; OU stationary distribution σ²_V = k_B T/m_p
    vector Vnew = VdragNew;
    if (brownianActive_)
    {
        const scalar varV = (kB*T_/mP_)*(1.0 - exp(-2.0*omega));
        const scalar sigV = sqrt(max(varV, scalar(0)));
        Vnew += sigV*xi;
    }

    return {Vnew, Vdisp};
}


// ************************************************************************* //
