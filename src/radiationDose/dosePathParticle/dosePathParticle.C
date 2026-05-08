/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dosePathParticle.H"
#include "dosePathCloud.H"
#include "constants.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(dosePathParticle, 0);
}
}


const Foam::NamedEnum<Foam::dose::dosePathParticle::endReason, 5>
Foam::dose::dosePathParticle::endReasonNames
{
    "active",
    "escaped",
    "timedOut",
    "stuck",
    "terminated"
};


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::dose::dosePathParticle::dosePathParticle
(
    const meshSearch& searchEngine,
    const vector& position,
    const label celli,
    label& nLocateBoundaryHits
)
:
    particle(searchEngine, position, celli, nLocateBoundaryHits),
    V_(Zero),
    V_disp_(Zero),
    D_(0),
    t_(0),
    endReason_(endReason::active),
    dispState_(),
    motionState_(),
    points_()
{}


Foam::dose::dosePathParticle::dosePathParticle
(
    Istream& is,
    bool readFields
)
:
    particle(is, readFields),
    V_(Zero),
    V_disp_(Zero),
    D_(0),
    t_(0),
    endReason_(endReason::active),
    dispState_(),
    motionState_(),
    points_()
{}


Foam::dose::dosePathParticle::dosePathParticle(const dosePathParticle& p)
:
    particle(p),
    V_(p.V_),
    V_disp_(p.V_disp_),
    D_(p.D_),
    t_(p.t_),
    endReason_(p.endReason_),
    dispState_(),
    motionState_(),
    points_(p.points_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dose::dosePathParticle::move
(
    dosePathCloud& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    // Already terminated by an earlier outer step? Keep the particle
    // in the cloud (so the function object can read its final state at
    // the end of execute()) but do no more tracking.
    if (endReason_ != endReason::active)
    {
        return true;
    }

    // Sample U at the current position; add a turbulent fluctuation.
    // The dispersion model owns its per-particle state (DRW eddy
    // lifetime + last fluctuation), which lives on this particle —
    // keeping each particle's stochastic trajectory independent of
    // every other particle's.
    const tetIndices tetIs = currentTetIndices(td.mesh);
    const vector Umean = td.UInterp().interpolate(coordinates(), tetIs);
    const vector uPrime =
        cloud.dispersion().fluctuation
        (
            *dispState_,
            position(td.mesh),
            cell(),
            cloud.dtMax(),
            td.rng()
        );

    // Pre-draw a unit-Gaussian Brownian sample if the motion model
    // wants one. The single draw is reused across the (potentially
    // CFL-shortened) re-advance below so Brownian remains a
    // single, well-defined random kick per outer step.
    vector xi(Zero);
    if (cloud.motion().needsBrownianSample())
    {
        const scalar twoPi = constant::mathematical::twoPi;
        const scalar r1 = max(td.rng().scalar01(), small);
        const scalar r2 = td.rng().scalar01();
        const scalar r3 = max(td.rng().scalar01(), small);
        const scalar r4 = td.rng().scalar01();
        const scalar mag1 = sqrt(-2.0*log(r1));
        const scalar mag3 = sqrt(-2.0*log(r3));
        xi = vector
        (
            mag1*cos(twoPi*r2),
            mag1*sin(twoPi*r2),
            mag3*cos(twoPi*r4)
        );
    }

    // Provisional advance over the full outer step to discover the
    // displacement velocity for the CFL bound. For tracer this is
    // exactly equivalent to the previous V_ = U + u' assignment.
    motionModel::stepResult r =
        cloud.motion().advance
        (
            *motionState_,
            V_,
            Umean,
            uPrime,
            cloud.dtMax(),
            xi
        );

    // CFL-bounded outer-step duration. cbrt(V[celli]) is a coarse
    // characteristic cell size; for the Sozzi case (mostly hex cells
    // ~1 mm) it correlates well with the actual face-to-face span.
    const scalar cellSize = cbrt(td.mesh.cellVolumes()[cell()]);
    const scalar VdispMag = mag(r.Vdisp);

    if (VdispMag < small)
    {
        endReason_ = endReason::stuck;
        return true;
    }

    const scalar dt =
        min(cloud.dtMax(), cloud.cflMax()*cellSize/VdispMag);

    // If CFL tightened the step, redo with the smaller dt so V (for
    // the next step) and V_disp (used right below) reflect the
    // shorter integration window. xi is intentionally re-used so the
    // Brownian draw maps to a single physical kick per outer step.
    if (dt < cloud.dtMax() - small)
    {
        r = cloud.motion().advance
            (*motionState_, V_, Umean, uPrime, dt, xi);
    }

    V_      = r.V;
    V_disp_ = r.Vdisp;

    // Reset stepFraction so the inner loop spans 0->1 of THIS outer
    // step (rather than the OF time step, which is irrelevant in the
    // function-object context).
    reset(0);

    // Inner loop: walk the particle through the dt budget, hitting
    // any number of cell faces / wall reflections along the way.
    // Termination of the loop happens when (a) stepFraction reaches 1
    // (the dt budget is consumed), (b) endReason_ becomes non-active
    // (escape patch, stuck, etc.), (c) the cloud asks us to stop
    // (keepParticle=false), or (d) a processor transfer is queued.
    while
    (
        stepFraction() < 1
     && endReason_ == endReason::active
     && td.keepParticle
     && td.sendToProc == -1
    )
    {
        const scalar sfrac = stepFraction();
        const scalar f = 1 - sfrac;

        const tetIndices tetIs_pre = currentTetIndices(td.mesh);
        const scalar G_pre =
            max(scalar(0), td.GInterp().interpolate(coordinates(), tetIs_pre));

        trackToAndHitFace(f*dt*V_disp_, f, cloud, td);

        const tetIndices tetIs_post = currentTetIndices(td.mesh);
        const scalar G_post =
            max(scalar(0), td.GInterp().interpolate(coordinates(), tetIs_post));

        const scalar actualDt = (stepFraction() - sfrac)*dt;
        D_ += 0.5*(G_pre + G_post)*actualDt*Wm2_s_to_mJcm2;
        t_ += actualDt;
    }

    // Soft termination checks. maxTime / maxDose of 0 disable.
    if (cloud.maxTime() > 0 && t_ >= cloud.maxTime())
    {
        endReason_ = endReason::timedOut;
    }
    else if (cloud.maxDose() > 0 && D_ >= cloud.maxDose())
    {
        endReason_ = endReason::terminated;
    }

    // Record the end-of-outer-step position. Skipped when the cloud
    // has storeTrack disabled — the only consumer of the stored
    // trajectory is the function-object's VTK writer, so when VTK
    // output is off there is no reason to grow points_ past its
    // initial seed entry. Memory then stays O(N_particles) instead
    // of O(N_particles × residence_time / dtMax), which matters for
    // long-trajectory production runs.
    if (cloud.storeTrack())
    {
        points_.append
        (
            trackPoint
            (
                position(td.mesh),
                t_,
                D_,
                cell()
            )
        );
    }

    return true;
}


void Foam::dose::dosePathParticle::hitWallPatch
(
    dosePathCloud& cloud,
    trackingData& td
)
{
    if (!cloud.wallReflection())
    {
        // No-reflection mode: any non-escape boundary hit terminates
        // the track. Useful for absorbing walls (e.g. UV reactor's
        // interior surfaces if treated as perfect absorbers).
        endReason_ = endReason::stuck;
        return;
    }

    // Specular reflection: V - 2*(V.n)*n. The barycentric coordinates
    // are unchanged (the particle stays exactly on the boundary face).
    // Reflect both V_disp_ (the inner-loop transport velocity, so the
    // remaining dt budget is consumed in the reflected direction) and
    // V_ (the true particle velocity, so the next outer step's
    // motion-model update starts from the post-reflection state).
    // For inertial particles the moment-of-hit V on the OU trajectory
    // is between V_old and V_; reflecting V_ is the tractable
    // approximation — within the same O(dt) family as the rest of
    // the inner-step truncation.
    const vector nw = normal(td.mesh);
    V_      -= 2.0*(V_      & nw)*nw;
    V_disp_ -= 2.0*(V_disp_ & nw)*nw;
}


void Foam::dose::dosePathParticle::hitBasicPatch
(
    dosePathCloud& cloud,
    trackingData& td
)
{
    const label patchi = patch(td.mesh);
    if (cloud.escapePatchIDs().found(patchi))
    {
        endReason_ = endReason::escaped;
    }
    else
    {
        // Some unknown patch type; treat as stuck. The base class's
        // hitBasicPatch sets keepParticle=false (delete on the spot),
        // which would lose the dose information — we keep it instead.
        endReason_ = endReason::stuck;
    }
}


void Foam::dose::dosePathParticle::readFields
(
    Foam::lagrangian::Cloud<dosePathParticle>& c
)
{
    // No-op: the dose tracker does not persist particles to disk.
    particle::readFields(c);
}


void Foam::dose::dosePathParticle::writeFields
(
    const Foam::lagrangian::Cloud<dosePathParticle>& c
)
{
    // No-op: the dose tracker does not persist particles to disk.
    // The function object writes its own CSV / summary directly.
    particle::writeFields(c);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::dose::operator<<
(
    Ostream& os,
    const dosePathParticle& p
)
{
    os  << static_cast<const particle&>(p) << token::SPACE
        << p.V_ << token::SPACE
        << p.V_disp_ << token::SPACE
        << p.D_ << token::SPACE
        << p.t_;
    return os;
}


// ************************************************************************* //
