/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "dosePathCloud.H"
#include "Pstream.H"

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(dosePathCloud, 0);
}
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::dose::dosePathCloud::dosePathCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    const scalar dtMax,
    const scalar cflMax,
    const labelHashSet& escapePatchIDs,
    const scalar maxTime,
    const scalar maxDose,
    const Switch wallReflection,
    const Switch storeTrack,
    autoPtr<dispersionModel> dispersion,
    autoPtr<motionModel> motion
)
:
    lagrangian::Cloud<dosePathParticle>(mesh, cloudName, false),
    mesh_(mesh),
    dtMax_(dtMax),
    cflMax_(cflMax),
    escapePatchIDs_(escapePatchIDs),
    maxTime_(maxTime),
    maxDose_(maxDose),
    wallReflection_(wallReflection),
    storeTrack_(storeTrack),
    dispersion_(std::move(dispersion)),
    motion_(std::move(motion))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dose::dosePathCloud::nActive() const
{
    label n = 0;
    forAllConstIter
    (
        typename Foam::lagrangian::Cloud<dosePathParticle>,
        *this,
        iter
    )
    {
        if (iter().active())
        {
            ++n;
        }
    }
    return n;
}


Foam::label Foam::dose::dosePathCloud::runToCompletion
(
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& GInterp,
    randomGenerator& parentRng,
    const label maxOuterSteps
)
{
    // In MPI mode the cloud iteration also has to flush per-rank
    // sendParticles[] queues across processor patches; that machinery
    // lives inside Cloud::move() and is not currently thread-safe, so
    // we keep MPI runs on the serial path and only OMP-parallelise
    // single-rank runs.
    const bool useOmp = !Pstream::parRun();

    // Serial fallback: one trackingData reused across outer steps
    // (preserves the rng stream so the result matches the pre-OMP
    // baseline for multi-rank or thread-disabled runs).
    dosePathParticle::trackingData serialTd
        (*this, UInterp, GInterp, parentRng);

    // OMP setup is hoisted out of the outer-step loop: per-thread
    // RNGs and trackingData persist for the lifetime of the run, so
    // each thread's RNG stream evolves naturally across outer steps
    // without redrawing seeds. The particle pointer vector is also
    // built once because the cloud is structurally stable — we
    // never delete particles (terminal endReasons are flagged
    // in-place) and single-rank runs have no processor handoff.
#ifdef _OPENMP
    const int nThreads = useOmp ? omp_get_max_threads() : 1;
#else
    const int nThreads = 1;
#endif

    std::vector<dosePathParticle*> particles;
    std::vector<randomGenerator> threadRngs;
    std::vector<dosePathParticle::trackingData> tds;

    if (useOmp)
    {
        particles.reserve(this->size());
        forAllIter
        (
            typename Foam::lagrangian::Cloud<dosePathParticle>,
            *this,
            iter
        )
        {
            particles.push_back(&iter());
        }

        // Derive per-thread seeds from independent draws against the
        // parent. Drawing once per thread (rather than seed = base + t)
        // is the canonical defensive pattern: even if the underlying
        // generator family has poor mixing for nearby seeds, the parent
        // stream's spacing between successive uniforms is the actual
        // seed gap. OF's Mersenne Twister has very strong seed mixing
        // so the "base + t" form was also fine in practice, but this
        // costs nothing extra and is the right pattern to copy elsewhere.
        threadRngs.reserve(nThreads);
        for (int t = 0; t < nThreads; ++t)
        {
            const label seed = label(parentRng.scalar01()*1.0e9);
            threadRngs.emplace_back(randomGenerator::seed(seed));
        }

        // trackingData has reference members; reserve+emplace_back
        // avoids reallocation that would invalidate copies.
        tds.reserve(nThreads);
        for (int t = 0; t < nThreads; ++t)
        {
            tds.emplace_back(*this, UInterp, GInterp, threadRngs[t]);
        }
    }

    label step = 0;
    while (step < maxOuterSteps)
    {
        // Global active-particle count: in MPI mode, ranks must enter
        // and exit the move() collective in lockstep. If only the
        // local count were checked, a rank that finished its particles
        // first would break out while peers still have live particles
        // to track, and the next move() collective on the lagging
        // ranks would hang.
        label active = nActive();
        if (Pstream::parRun())
        {
            reduce(active, sumOp<label>());
        }
        if (active == 0)
        {
            break;
        }
        if (debug && (step % 100 == 0))
        {
            Info<< "  outer step " << step
                << ": " << active << " active" << endl;
        }
        if (useOmp)
        {
            moveOmpStep(particles, tds);
        }
        else
        {
            Cloud<dosePathParticle>::move(*this, serialTd);
        }
        ++step;
    }
    return step;
}


void Foam::dose::dosePathCloud::moveOmpStep
(
    std::vector<dosePathParticle*>& particles,
    std::vector<dosePathParticle::trackingData>& tds
)
{
    if (particles.empty()) return;

    // Static schedule with default chunking gives each thread a
    // contiguous block of indices: the result is a deterministic
    // function of (parentSeed, nThreads, particle order), and thread
    // T always processes the same slice every outer step. The
    // dispersion model state on each particle is self-contained;
    // the only cross-particle shared state is the RNG, which is
    // per-thread, so there is no coordination inside the loop.
    const long n = static_cast<long>(particles.size());

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (long i = 0; i < n; ++i)
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        particles[i]->move(*this, tds[tid]);
    }
}


// ************************************************************************* //
