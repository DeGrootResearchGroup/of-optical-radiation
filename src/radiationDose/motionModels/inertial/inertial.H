/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::inertial

Description
    Inertial point-particle motion model. Integrates the linear-drag
    Langevin equation

        dV/dt = (U_seen - V)/τ_p + a_g + a_B(t)

    over one outer step using the Ornstein-Uhlenbeck exact update
    (analytical for piecewise-constant U_seen and a_g over Δt; gives
    the correct stationary variance for the Brownian kick).

    Symbols and equations
    ---------------------
    U_seen = U + u'                                    (m/s)
    a_g    = (1 - ρ_f/ρ_p) g                            (m/s²)   (gravity off ⇒ 0)
    τ_p    = drag-model response time                  (s)
    V_eq   = U_seen + a_g·τ_p                          (m/s)
    ω      = Δt / τ_p
    φ      = (1 - exp(-ω))/ω    (→ 1 as ω→0, → 0 as ω→∞)

    End-of-step velocity (drag only):
        V_drag(Δt) = V_eq + (V₀ - V_eq) · exp(-ω)

    Mean velocity over [0, Δt] (Δx/Δt; used as V_disp by the caller
    for trackToAndHitFace):
        V_disp     = V_eq + (V₀ - V_eq) · φ

    Brownian kick added to V at end of step (FDT-tied to the same τ_p):
        ΔV_B = √(k_B·T/m_p · (1 - exp(-2ω))) · ξ,   ξ ~ N(0, I)

    Per-component sigma is the OU stationary √variance scaled by the
    relaxation factor; recovers the Maxwell-Boltzmann stationary
    distribution V ~ N(V_eq, k_B·T/m_p · I) as ω → ∞.

    The position-noise contribution from the velocity Brownian motion
    over [0, Δt] is dropped from V_disp. This under-counts mean-square
    displacement on time-scales Δt ≲ τ_p (where the missing term is
    O(τ_p²) per step) but recovers the correct long-time diffusion
    coefficient D = k_B·T·τ_p/m_p as Brownian-driven V correlations
    decay between steps. Add the σ_x term to V_disp if a driver case
    needs sub-τ_p Brownian displacement statistics.

    Schiller-Naumann coupling
    -------------------------
    For Re_p-corrected drag, τ_p depends on |U_seen - V|. We evaluate
    τ_p once per outer step at the start-of-step V; this is the linear
    drag coefficient that the OU update treats as constant over Δt.
    Lagging is consistent with the OU linearisation and avoids any
    fixed-point iteration.

    Configuration
    -------------
        motion
        {
            type     inertial;
            rhoP     1050;     // kg/m^3
            dP       50e-6;    // m
            rhoF     1000;     // kg/m^3 (default 1000, water)
            muF      1e-3;     // Pa·s   (default 1e-3, water at 293 K)

            drag
            {
                type    schillerNaumann;   // or stokes
            }

            gravity
            {
                active  true;
                value   (0 -9.81 0);
            }

            brownian
            {
                active  false;
                T       293.15;             // K
            }
        }

\*---------------------------------------------------------------------------*/

#ifndef dose_inertial_H
#define dose_inertial_H

#include "Switch.H"
#include "autoPtr.H"
#include "dragModel.H"
#include "motionModel.H"

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class inertial Declaration
\*---------------------------------------------------------------------------*/

class inertial
:
    public motionModel
{
    // Private Data

        autoPtr<dragModel> drag_;

        scalar rhoP_;          // kg/m^3
        scalar dP_;            // m
        scalar rhoF_;          // kg/m^3
        scalar muF_;           // Pa·s
        scalar nuF_;           // m^2/s = muF_/rhoF_
        scalar mP_;            // kg, (pi/6)·rhoP·dP^3

        Switch  gravityActive_;
        vector  gravity_;      // m/s^2 (signed)

        Switch  brownianActive_;
        scalar  T_;            // K


public:

    //- Runtime type information
    TypeName("inertial");


    // Constructors

        inertial(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~inertial() = default;


    // Member Functions

        virtual bool needsBrownianSample() const
        {
            return brownianActive_;
        }

        virtual stepResult advance
        (
            State& state,
            const vector& V,
            const vector& Umean,
            const vector& uPrime,
            scalar dt,
            const vector& xi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
