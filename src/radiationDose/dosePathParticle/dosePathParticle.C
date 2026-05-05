/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dosePathParticle.H"
#include "dosePathCloud.H"
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
    D_(0),
    t_(0),
    endReason_(endReason::active),
    dispState_(),
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
    D_(0),
    t_(0),
    endReason_(endReason::active),
    dispState_(),
    points_()
{}


Foam::dose::dosePathParticle::dosePathParticle(const dosePathParticle& p)
:
    particle(p),
    V_(p.V_),
    D_(p.D_),
    t_(p.t_),
    endReason_(p.endReason_),
    dispState_(),
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
    V_ = Umean + uPrime;

    if (mag(V_) < small)
    {
        endReason_ = endReason::stuck;
        return true;
    }

    // CFL-bounded outer-step duration. cbrt(V[celli]) is a coarse
    // characteristic cell size; for the Sozzi case (mostly hex cells
    // ~1 mm) it correlates well with the actual face-to-face span.
    const scalar cellSize = cbrt(td.mesh.cellVolumes()[cell()]);
    const scalar dt = min(cloud.dtMax(), cloud.cflMax()*cellSize/mag(V_));

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

        trackToAndHitFace(f*dt*V_, f, cloud, td);

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

    // Record the end-of-outer-step position. The trajectory is
    // currently always stored; a future storeFullTrack switch can
    // make this conditional.
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
    // are unchanged (the particle stays exactly on the boundary face)
    // and the inner loop continues with the new V_ for the remainder
    // of the dt budget.
    const vector nw = normal(td.mesh);
    V_ -= 2.0*(V_ & nw)*nw;
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
        << p.D_ << token::SPACE
        << p.t_;
    return os;
}


// ************************************************************************* //
