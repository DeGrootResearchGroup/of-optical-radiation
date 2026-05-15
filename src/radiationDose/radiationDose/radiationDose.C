/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "radiationDose.H"
#include "IOobject.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "meshSearch.H"
#include "polyBoundaryMesh.H"
#include "volFields.H"

#include <cmath>
#include <cstdio>
#include <limits>

namespace
{
    // Clamp values whose magnitude falls below float32's smallest
    // normal (~1.18e-38) to exactly 0 before VTK emission. ParaView's
    // legacy-ASCII reader loses sync with the declared array size when
    // it encounters subnormal floats -- empirically the parser bails
    // mid-array on the next-array `SCALARS` keyword and the file
    // becomes unreadable past the first one or two arrays. Dose / time
    // accumulators routinely produce O(1e-20..1e-45) values in cells
    // with negligible G (e.g. shadowed seed positions), and those
    // values would have rounded to 0 once stored as the SCALARS-
    // declared `float` type anyway -- so this flush-to-zero is
    // information-preserving for the VTK file format.
    inline Foam::scalar vtkSafeFloat(const Foam::scalar x)
    {
        return
            std::abs(x) < std::numeric_limits<float>::min()
          ? Foam::scalar(0)
          : x;
    }

    // Zero-padded numeric suffix for per-batch filenames. Five
    // digits is enough for 99999 batches; at the dictionary default
    // batchSize of 0 (disabled) the limit doesn't bite.
    inline Foam::word padIndex(const Foam::label i, const int width = 5)
    {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%0*ld", width, long(i));
        return Foam::word(buf);
    }
}

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(radiationDose, 0);
    addToRunTimeSelectionTable(functionObject, radiationDose, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::radiationDose::radiationDose
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mode_("steady"),
    UName_("U"),
    GName_("G"),
    randomSeed_(0),
    maxTime_(0),
    maxDose_(0),
    wallReflection_(true),
    dtMax_(0.01),
    cflMax_(0.5),
    maxOuterSteps_(100000),
    writeVtk_(true),
    trajectoryStride_(1),
    vtkMinDose_(-1),
    vtkMaxDose_(-1),
    batchSize_(0),
    seeded_(false),
    emitted_(false)
{
    read(dict);
}


Foam::functionObjects::radiationDose::~radiationDose()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::radiationDose::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    mode_ = dict.lookupOrDefault<word>("mode", word("steady"));
    if (mode_ != "steady" && mode_ != "unsteady")
    {
        FatalErrorInFunction
            << "Unknown mode \"" << mode_ << "\". Expected"
            << " \"steady\" or \"unsteady\"."
            << exit(FatalError);
    }

    UName_ = dict.lookupOrDefault<word>("U", "U");
    GName_ = dict.lookupOrDefault<word>("fluenceRate", "G");
    randomSeed_ = dict.lookupOrDefault<label>("seed", 0);

    seeding_ = dose::seedingModel::New(dict.subDict("seeding"), mesh_);
    dispersionDict_ = dict.subDict("dispersion");

    if (dict.found("motion"))
    {
        motionDict_ = dict.subDict("motion");
    }
    else
    {
        // Default: fluid-tracer motion (V = U + u'). Matches the v0.3
        // behaviour and keeps existing tutorials minimal.
        motionDict_.clear();
        motionDict_.add("type", word("tracer"));
    }

    const dictionary& termDict = dict.subDict("termination");
    const wordList escapePatches(termDict.lookup("escapePatches"));
    escapePatchIDs_.clear();
    forAll(escapePatches, i)
    {
        // Each rank looks the patch up locally. For mesh boundary
        // patches (inlet/outlet/walls/etc.) the index is identical on
        // every rank because OF mirrors the full boundary list on each
        // decomposed rank (zero-face stubs where the patch isn't
        // owned). A missing name on a rank therefore means a real
        // typo, not a decomposition artefact -- but in parallel we
        // reduce-then-fail so all ranks abort with the same message
        // instead of just the first one to notice.
        const label patchi =
            mesh_.boundaryMesh().findIndex(escapePatches[i]);

        label maxPatchi = patchi;
        if (Pstream::parRun())
        {
            reduce(maxPatchi, maxOp<label>());
        }
        if (maxPatchi < 0)
        {
            FatalErrorInFunction
                << "escapePatch " << escapePatches[i]
                << " not found on any rank"
                << exit(FatalError);
        }
        if (patchi >= 0)
        {
            escapePatchIDs_.insert(patchi);
        }
    }
    maxTime_ = termDict.lookupOrDefault<scalar>("maxTime", 0);
    maxDose_ = termDict.lookupOrDefault<scalar>("maxDose", 0);
    wallReflection_ =
        termDict.lookupOrDefault<Switch>("wallReflection", true);

    const dictionary& intDict = dict.subDict("integration");
    dtMax_ = readScalar(intDict.lookup("dtMax"));
    cflMax_ = readScalar(intDict.lookup("cflMax"));
    // Old name "maxStepsPerTrack" was a per-track inner-loop cap; the
    // new tracker has only an outer-step cap (number of dtMax steps).
    // Accept both names for now to keep tutorials working.
    maxOuterSteps_ =
        intDict.lookupOrDefault<label>
        (
            "maxOuterSteps",
            intDict.lookupOrDefault<label>("maxStepsPerTrack", 100000)
        );

    const dictionary& outDict = dict.subDict("output");
    kInact_ = outDict.lookupOrDefault<scalarList>("kInact", scalarList());
    forAll(kInact_, i)
    {
        const scalar k = kInact_[i];
        if (k <= 0)
        {
            FatalErrorInFunction
                << "kInact[" << i << "] = " << k
                << " cm^2/mJ is not positive. Each entry must be > 0;"
                << " the dose-response model exp(-k*D) requires it."
                << exit(FatalError);
        }
        // Soft band on typical microbial UV inactivation rate constants
        // (cm^2/mJ): ~5e-3 for spores up to ~1.5 for highly UV-sensitive
        // bacteria. A value far outside this band almost certainly
        // signals a unit mix-up (e.g. m^2/J = 10 cm^2/mJ).
        if (k < 1e-3 || k > 100)
        {
            WarningInFunction
                << "kInact[" << i << "] = " << k
                << " cm^2/mJ is well outside the typical microbial"
                << " range (~5e-3 to 1.5 cm^2/mJ). Check the units."
                << endl;
        }
    }
    writeVtk_ = outDict.lookupOrDefault<Switch>("writeVtk", true);

    trajectoryStride_ =
        outDict.lookupOrDefault<label>("trajectoryStride", 1);
    if (trajectoryStride_ < 1)
    {
        FatalErrorInFunction
            << "trajectoryStride must be >= 1 (got "
            << trajectoryStride_ << ")."
            << exit(FatalError);
    }

    vtkMinDose_ = outDict.lookupOrDefault<scalar>("vtkMinDose", -1);
    vtkMaxDose_ = outDict.lookupOrDefault<scalar>("vtkMaxDose", -1);
    if (vtkMinDose_ >= 0 && vtkMaxDose_ >= 0 && vtkMaxDose_ < vtkMinDose_)
    {
        FatalErrorInFunction
            << "vtkMaxDose (" << vtkMaxDose_ << ") < vtkMinDose ("
            << vtkMinDose_ << "). Either bound can be disabled by"
            << " setting it to a negative value."
            << exit(FatalError);
    }

    batchSize_ = outDict.lookupOrDefault<label>("batchSize", 0);
    if (batchSize_ < 0)
    {
        FatalErrorInFunction
            << "batchSize must be >= 0 (got " << batchSize_
            << "). Use 0 to disable batching."
            << exit(FatalError);
    }
    if (batchSize_ > 0 && Pstream::parRun())
    {
        // Batching slices the global seed list into chunks, which
        // doesn't compose with the per-rank seed distribution that
        // patchInjection / pointInjection produce in parallel.
        // Single-rank cases (incl. multi-thread OMP) work fine.
        FatalErrorInFunction
            << "batchSize > 0 is not supported in parallel runs."
            << " Run single-rank (OMP threading via OMP_NUM_THREADS"
            << " still works) or disable batching."
            << exit(FatalError);
    }
    if (batchSize_ > 0 && mode_ == "unsteady")
    {
        // Batching is a steady-mode memory-bound device: each batch
        // is a fresh cloud run to completion. Unsteady mode has a
        // single cohort spanning the host time loop, so batching
        // does not apply.
        FatalErrorInFunction
            << "batchSize > 0 is not supported in unsteady mode."
            << " Set batchSize 0 (the default) or run in steady mode."
            << exit(FatalError);
    }

    return true;
}


Foam::wordList Foam::functionObjects::radiationDose::fields() const
{
    DynamicList<word> f;
    f.append(UName_);
    f.append(GName_);

    // Static dispatch: each concrete dispersion / motion model
    // registers a wordList(*)(const dictionary&) at static-init time,
    // and the base class's requiredFields(dict) looks it up by type
    // name. No throwaway model is constructed -- in particular, the
    // inertial motion model doesn't need to build its drag-model
    // sub-tree just to tell us it needs no extra fields.
    const wordList dispExtra =
        dose::dispersionModel::requiredFields(dispersionDict_);
    forAll(dispExtra, i) { f.append(dispExtra[i]); }

    const wordList motExtra =
        dose::motionModel::requiredFields(motionDict_);
    forAll(motExtra, i) { f.append(motExtra[i]); }

    return wordList(f);
}


// Per-rank gather payload. Per-track scalars are kept as a List per
// rank; per-vertex data is flattened (concatenated across tracks) and
// sliced via the matching nVerts list. nVerts[r][i] is 0 for tracks
// that should not appear in the VTK file (writeVtk_ off, or fewer
// than 2 vertices); the CSV/summary writers ignore it.
struct Foam::functionObjects::radiationDose::GatheredTracks
{
    // Per-track summary, indexed [proc][trackOnProc]
    List<List<label>>  origId;
    List<List<label>>  origProc;
    List<List<label>>  endReason;
    List<List<scalar>> tEnd;
    List<List<scalar>> dose;
    List<List<vector>> xEnd;

    // Per-track vertex count (0 if the track is excluded from VTK)
    List<List<label>>  nVerts;

    // Per-vertex data, concatenated across tracks within each rank
    List<List<vector>> vertX;
    List<List<scalar>> vertT;
    List<List<scalar>> vertD;
    List<List<label>>  vertC;

    // Globally-unique trackId per track, indexed [proc][trackOnProc].
    // Computed by runBatch() right after gather: equal to the
    // per-particle origId offset by the running aggregator size at
    // the start of this batch. Used as the VTK CELL_DATA trackId
    // scalar and propagated into the AggregatedTracks for the CSV.
    List<List<label>>  globalTrackId;
};


// Cross-batch accumulator. Flat 1-D, master-only, summary fields plus
// per-batch VTK filenames. Populated incrementally inside execute()
// via appendBatch() as each batch finishes; consumed by write() for
// the aggregated CSV + summary + PVD wrapper.
struct Foam::functionObjects::radiationDose::AggregatedTracks
{
    DynamicList<label>  trackId;     // global seed index
    DynamicList<label>  endReason;   // integer enum cast
    DynamicList<scalar> tEnd;
    DynamicList<scalar> dose;
    DynamicList<vector> xEnd;

    DynamicList<word>   vtkFile;     // basename of each batch's .vtk
};


bool Foam::functionObjects::radiationDose::execute()
{
    // The host driver registers U automatically (every fluid solver
    // owns the velocity field), but `G` is a radiation field that
    // most flow solvers don't know about. foamPostProcess auto-loads
    // any volScalarField it finds in the time directory, but foamRun
    // / foamMultiRun load only the fields their solver declares. To
    // keep the FO self-contained in either driver, fall back to a
    // disk read via findInstance() the first time G isn't in the
    // registry. The new VolField self-registers with the mesh
    // registry and persists for the rest of the run.
    if (!mesh_.foundObject<volScalarField>(GName_))
    {
        const fileName instance =
            time_.findInstance(mesh_.dbDir(), GName_);
        Info<< type() << ": loading " << GName_
            << " from " << instance << endl;
        new volScalarField
        (
            IOobject
            (
                GName_,
                instance,
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
    }

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& G =
        mesh_.lookupObject<volScalarField>(GName_);

    if (mode_ == "unsteady")
    {
        return executeUnsteady(U, G);
    }

    randomGenerator rng{randomGenerator::seed(randomSeed_)};

    // Seed all particles up front. The seedingModel returns this
    // rank's share; in parallel runs that's only the seeds whose
    // injection face / cell is owned locally, so nLocal varies across
    // ranks while nGlobal is the sum every rank agrees on.
    Info<< type() << ": seeding particles..." << endl;
    const List<dose::trackPoint> allSeeds = seeding_->seed(rng);
    const label nLocal = allSeeds.size();
    label nGlobal = nLocal;
    if (Pstream::parRun())
    {
        reduce(nGlobal, sumOp<label>());
    }
    Info<< type() << ": seeded " << nGlobal << " particles" << endl;

    // Decide batching layout. batchSize == 0, or batchSize >= nGlobal,
    // collapses to a single batch keeping the original filename
    // (trajectories.vtk) and skipping the PVD wrapper. The loop count
    // is set from the *global* total because cloud construction inside
    // runBatch is a collective and every rank must enter every batch
    // in lockstep -- a rank whose local seed share happens to be empty
    // still has to participate in the batch's MPI alltoall.
    const label effBatch =
        (batchSize_ <= 0 || batchSize_ >= nGlobal) ? nGlobal : batchSize_;
    const label nBatches =
        (nGlobal > 0) ? (nGlobal + effBatch - 1) / effBatch : 0;
    const bool batched = (nBatches > 1);

    fileName outDir;
    if (Pstream::master())
    {
        outDir =
            time_.globalPath()/"postProcessing"/name()/time_.name();
        mkDir(outDir);
        aggregated_.reset(new AggregatedTracks);
    }

    if (nBatches == 0)
    {
        // Edge case: no seeds anywhere. Empty CSV / summary still get
        // written from write() against the (empty) aggregator.
        return true;
    }

    for (label b = 0; b < nBatches; ++b)
    {
        const word vtkFile = batched
            ? word("trajectories_" + padIndex(b + 1) + ".vtk")
            : word("trajectories.vtk");

        // Two slicing regimes:
        //   * Batched -- single-rank only (read() enforces it). The
        //     `allSeeds` list is the full global seed array; we slice
        //     it by global index.
        //   * Unbatched -- nBatches == 1. Each rank passes its full
        //     local seed list straight through, so the union across
        //     ranks reconstitutes the global population. This is the
        //     only path that survives a multi-rank run.
        List<dose::trackPoint> seeds;
        if (batched)
        {
            const label start = b*effBatch;
            const label end   = min(start + effBatch, nLocal);
            const label n     = end - start;
            Info<< type() << ": batch " << (b + 1) << "/" << nBatches
                << " (" << n << " particles, global indices "
                << start << ".." << (end - 1) << ")" << endl;
            seeds.setSize(n);
            for (label i = 0; i < n; ++i)
            {
                seeds[i] = allSeeds[start + i];
            }
        }
        else
        {
            seeds = allSeeds;
        }

        runBatch(seeds, vtkFile, outDir, U, G, rng);
    }

    return true;
}


void Foam::functionObjects::radiationDose::runBatch
(
    const List<dose::trackPoint>& seeds,
    const word& vtkFileName,
    const fileName& outDir,
    const volVectorField& U,
    const volScalarField& G,
    randomGenerator& rng
)
{
    // Build the dispersion and motion models fresh per batch so any
    // internal counters / RNG state start clean.
    autoPtr<dose::dispersionModel> dispersion =
        dose::dispersionModel::New(dispersionDict_, mesh_);
    autoPtr<dose::motionModel> motion =
        dose::motionModel::New(motionDict_, mesh_);

    cloud_.reset
    (
        new dose::dosePathCloud
        (
            mesh_,
            "doseCloud",
            dtMax_,
            cflMax_,
            escapePatchIDs_,
            maxTime_,
            maxDose_,
            wallReflection_,
            writeVtk_,
            trajectoryStride_,
            std::move(dispersion),
            std::move(motion)
        )
    );

    const meshSearch& ms = meshSearch::New(mesh_);
    label nLocateBoundaryHits = 0;
    forAll(seeds, i)
    {
        dose::dosePathParticle* p =
            new dose::dosePathParticle
            (
                ms,
                seeds[i].xd(),
                seeds[i].celli,
                nLocateBoundaryHits
            );
        p->setDispState(cloud_->dispersion().newState());
        p->setMotionState(cloud_->motion().newState());
        if (writeVtk_)
        {
            p->appendPoint(seeds[i]);
        }
        cloud_->addParticle(p);
    }

    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> GInterp(G);

    Info<< type() << ": integrating..." << endl;
    const label nSteps =
        cloud_->runToCompletion(UInterp, GInterp, rng, maxOuterSteps_);

    label nEscaped = 0, nTimedOut = 0, nStuck = 0, nTerm = 0, nActive = 0;
    forAllConstIter
    (
        typename lagrangian::Cloud<dose::dosePathParticle>,
        *cloud_,
        iter
    )
    {
        switch (iter().end())
        {
            case dose::dosePathParticle::endReason::escaped:    ++nEscaped;  break;
            case dose::dosePathParticle::endReason::timedOut:   ++nTimedOut; break;
            case dose::dosePathParticle::endReason::stuck:      ++nStuck;    break;
            case dose::dosePathParticle::endReason::terminated: ++nTerm;     break;
            case dose::dosePathParticle::endReason::active:     ++nActive;   break;
        }
    }
    Info<< type() << ": done in " << nSteps << " outer steps. "
        << nEscaped  << " escaped, "
        << nTimedOut << " timed out, "
        << nStuck    << " stuck, "
        << nTerm     << " hit maxDose, "
        << nActive   << " still active (capped by maxOuterSteps)"
        << endl;

    GatheredTracks g = gatherTracks();

    if (Pstream::master())
    {
        // Stamp the globally-unique trackId onto each track. The
        // per-particle origId resets to 0 with every new cloud, so we
        // recover a globally consistent index by offsetting by the
        // running aggregator size before appending this batch.
        const label trackOffset = aggregated_->trackId.size();
        forAll(g.globalTrackId, r)
        {
            const label nr = g.origId[r].size();
            g.globalTrackId[r].setSize(nr);
            for (label i = 0; i < nr; ++i)
            {
                g.globalTrackId[r][i] = trackOffset + g.origId[r][i];
            }
        }

        if (writeVtk_)
        {
            writeVtkTrajectories(g, outDir/vtkFileName);
        }
        appendBatch(g, vtkFileName);
    }

    // Drop the cloud before the next batch so its per-particle
    // trajectory storage is released (the gather already pulled
    // everything needed for CSV/summary/VTK to rank 0).
    cloud_.clear();
}


void Foam::functionObjects::radiationDose::appendBatch
(
    const GatheredTracks& g,
    const word& vtkFileName
)
{
    AggregatedTracks& a = aggregated_();
    forAll(g.origId, r)
    {
        const label nr = g.origId[r].size();
        for (label i = 0; i < nr; ++i)
        {
            a.trackId.append(g.globalTrackId[r][i]);
            a.endReason.append(g.endReason[r][i]);
            a.tEnd.append(g.tEnd[r][i]);
            a.dose.append(g.dose[r][i]);
            a.xEnd.append(g.xEnd[r][i]);
        }
    }
    if (writeVtk_)
    {
        a.vtkFile.append(vtkFileName);
    }
}


bool Foam::functionObjects::radiationDose::executeUnsteady
(
    const volVectorField& U,
    const volScalarField& G
)
{
    if (!seeded_)
    {
        // Restart guard: in unsteady mode the cohort lives only in
        // memory, so a run that starts mid-simulation has no way to
        // pick up where the previous run left off. Fail loud rather
        // than silently re-seeding a "fresh" cohort at the restart
        // time -- the user almost certainly does not want that.
        //
        // The first execute() fires AFTER the first time step
        // completes (default OF function-object semantics; the FO
        // does not opt in to executeAtStart), so tNow at first
        // call is roughly tStart + deltaT for a fresh run. Allow
        // up to 2 deltaT of slack to absorb both that offset and
        // any cumulative-write-interval cadence -- a real restart
        // from a checkpoint at t >> tStart blows past this margin
        // by orders of magnitude.
        const scalar tNow = time_.value();
        const scalar tStart = time_.startTime().value();
        const scalar dt = time_.deltaT().value();
        if (tNow > tStart + 2*dt + small)
        {
            FatalErrorInFunction
                << "unsteady mode does not support restart. The run"
                << " is at t = " << tNow << " s but the configured"
                << " startTime is " << tStart << " s. Restart with"
                << " startTime equal to the configured startTime,"
                << " or use steady mode against the latest stored"
                << " time."
                << exit(FatalError);
        }

        rng_.reset
        (
            new randomGenerator
            (
                randomGenerator::seed(randomSeed_)
            )
        );

        Info<< type() << ": seeding particles (unsteady mode)..." << endl;
        const List<dose::trackPoint> seeds = seeding_->seed(*rng_);
        // Same nLocal/nGlobal accounting as steady mode: seeds is the
        // local share; report the cross-rank sum so the log lines up
        // with the actual cohort size in parallel runs.
        label nGlobal = seeds.size();
        if (Pstream::parRun())
        {
            reduce(nGlobal, sumOp<label>());
        }
        Info<< type() << ": seeded " << nGlobal << " particles" << endl;

        if (Pstream::master())
        {
            aggregated_.reset(new AggregatedTracks);
        }

        autoPtr<dose::dispersionModel> dispersion =
            dose::dispersionModel::New(dispersionDict_, mesh_);
        autoPtr<dose::motionModel> motion =
            dose::motionModel::New(motionDict_, mesh_);

        cloud_.reset
        (
            new dose::dosePathCloud
            (
                mesh_,
                "doseCloud",
                dtMax_,
                cflMax_,
                escapePatchIDs_,
                maxTime_,
                maxDose_,
                wallReflection_,
                writeVtk_,
                trajectoryStride_,
                std::move(dispersion),
                std::move(motion)
            )
        );

        const meshSearch& ms = meshSearch::New(mesh_);
        label nLocateBoundaryHits = 0;
        forAll(seeds, i)
        {
            dose::dosePathParticle* p =
                new dose::dosePathParticle
                (
                    ms,
                    seeds[i].xd(),
                    seeds[i].celli,
                    nLocateBoundaryHits
                );
            p->setDispState(cloud_->dispersion().newState());
            p->setMotionState(cloud_->motion().newState());
            if (writeVtk_)
            {
                p->appendPoint(seeds[i]);
            }
            cloud_->addParticle(p);
        }

        seeded_ = true;
    }

    // Build interpolators fresh each call: U and G may have been
    // updated by the host solver since the last execute().
    interpolationCellPoint<vector> UInterp(U);
    interpolationCellPoint<scalar> GInterp(G);

    const scalar dt = time_.deltaT().value();
    if (dt > 0)
    {
        cloud_->runForDuration
        (
            dt,
            UInterp,
            GInterp,
            *rng_,
            maxOuterSteps_
        );
    }

    return true;
}


void Foam::functionObjects::radiationDose::flushUnsteady()
{
    if (emitted_ || !cloud_.valid())
    {
        return;
    }

    label nEscaped = 0, nTimedOut = 0, nStuck = 0, nTerm = 0, nActive = 0;
    forAllConstIter
    (
        typename lagrangian::Cloud<dose::dosePathParticle>,
        *cloud_,
        iter
    )
    {
        switch (iter().end())
        {
            case dose::dosePathParticle::endReason::escaped:    ++nEscaped;  break;
            case dose::dosePathParticle::endReason::timedOut:   ++nTimedOut; break;
            case dose::dosePathParticle::endReason::stuck:      ++nStuck;    break;
            case dose::dosePathParticle::endReason::terminated: ++nTerm;     break;
            case dose::dosePathParticle::endReason::active:     ++nActive;   break;
        }
    }
    Info<< type() << ": unsteady run ending. "
        << nEscaped  << " escaped, "
        << nTimedOut << " timed out, "
        << nStuck    << " stuck, "
        << nTerm     << " hit maxDose, "
        << nActive   << " still active at simulation end"
        << endl;

    GatheredTracks g = gatherTracks();

    if (Pstream::master())
    {
        const label trackOffset = aggregated_->trackId.size();
        forAll(g.globalTrackId, r)
        {
            const label nr = g.origId[r].size();
            g.globalTrackId[r].setSize(nr);
            for (label i = 0; i < nr; ++i)
            {
                g.globalTrackId[r][i] = trackOffset + g.origId[r][i];
            }
        }

        const fileName outDir =
            time_.globalPath()/"postProcessing"/name()/time_.name();
        mkDir(outDir);

        if (writeVtk_)
        {
            writeVtkTrajectories(g, outDir/"trajectories.vtk");
        }
        appendBatch(g, "trajectories.vtk");

        const AggregatedTracks& a = aggregated_();
        writeDoseCsv(a, outDir);
        writeSummary(a, outDir);
    }

    emitted_ = true;
    cloud_.clear();
}


bool Foam::functionObjects::radiationDose::write()
{
    // Unsteady mode defers all output to end-of-run; write() called
    // at intermediate write intervals is a no-op (intermediate state
    // is meaningless for a single in-flight cohort).
    if (mode_ == "unsteady")
    {
        return true;
    }

    if (!aggregated_.valid())
    {
        return true;
    }

    if (Pstream::master())
    {
        // globalPath() strips the processorN suffix in parallel runs;
        // in serial it's identical to path(). aggregated_ is already
        // global (rank 0 holds all batches' data after each gather),
        // so the CSV / summary / PVD belong in the case-root
        // postProcessing tree, not under any processor directory.
        const fileName outDir =
            time_.globalPath()/"postProcessing"/name()/time_.name();
        mkDir(outDir);

        const AggregatedTracks& a = aggregated_();
        writeDoseCsv(a, outDir);
        writeSummary(a, outDir);

        // PVD wrapper only when batching produced more than one VTK
        // file; a single-batch run keeps the original single-file
        // layout for backward compatibility.
        if (writeVtk_ && a.vtkFile.size() > 1)
        {
            writePvdWrapper(a, outDir);
        }
    }
    return true;
}


bool Foam::functionObjects::radiationDose::end()
{
    // Unsteady mode flushes at end-of-run. Steady mode is fully
    // driven by execute()/write() pairs from foamPostProcess so
    // end() has nothing to do.
    if (mode_ == "unsteady")
    {
        flushUnsteady();
    }
    return fvMeshFunctionObject::end();
}


Foam::functionObjects::radiationDose::GatheredTracks
Foam::functionObjects::radiationDose::gatherTracks() const
{
    GatheredTracks g;
    const label nProcs = Pstream::nProcs();
    g.origId.setSize(nProcs);
    g.origProc.setSize(nProcs);
    g.endReason.setSize(nProcs);
    g.tEnd.setSize(nProcs);
    g.dose.setSize(nProcs);
    g.xEnd.setSize(nProcs);
    g.nVerts.setSize(nProcs);
    g.vertX.setSize(nProcs);
    g.vertT.setSize(nProcs);
    g.vertD.setSize(nProcs);
    g.vertC.setSize(nProcs);
    g.globalTrackId.setSize(nProcs);

    const label me = Pstream::myProcNo();
    const label nLocal = cloud_->size();
    g.origId[me].setSize(nLocal);
    g.origProc[me].setSize(nLocal);
    g.endReason[me].setSize(nLocal);
    g.tEnd[me].setSize(nLocal);
    g.dose[me].setSize(nLocal);
    g.xEnd[me].setSize(nLocal);
    g.nVerts[me].setSize(nLocal);

    DynamicList<vector> vx;
    DynamicList<scalar> vt;
    DynamicList<scalar> vd;
    DynamicList<label>  vc;

    label k = 0;
    forAllConstIter
    (
        typename lagrangian::Cloud<dose::dosePathParticle>,
        *cloud_,
        iter
    )
    {
        const dose::dosePathParticle& p = iter();
        g.origId[me][k]    = p.origId();
        g.origProc[me][k]  = p.origProc();
        g.endReason[me][k] = label(p.end());
        g.tEnd[me][k]      = p.t();
        g.dose[me][k]      = p.D();
        g.xEnd[me][k]      = p.position(mesh_);

        const DynamicList<dose::trackPoint>& pts = p.points();
        const bool doseInRange =
            (vtkMinDose_ < 0 || p.D() >= vtkMinDose_)
         && (vtkMaxDose_ < 0 || p.D() <= vtkMaxDose_);
        const label nv =
            (writeVtk_ && pts.size() >= 2 && doseInRange) ? pts.size() : 0;
        g.nVerts[me][k] = nv;
        for (label j = 0; j < nv; ++j)
        {
            vx.append(pts[j].xd());
            vt.append(pts[j].t);
            vd.append(pts[j].D);
            vc.append(pts[j].celli);
        }
        ++k;
    }

    g.vertX[me] = vx;
    g.vertT[me] = vt;
    g.vertD[me] = vd;
    g.vertC[me] = vc;

    Pstream::gatherList(g.origId);
    Pstream::gatherList(g.origProc);
    Pstream::gatherList(g.endReason);
    Pstream::gatherList(g.tEnd);
    Pstream::gatherList(g.dose);
    Pstream::gatherList(g.xEnd);
    Pstream::gatherList(g.nVerts);
    Pstream::gatherList(g.vertX);
    Pstream::gatherList(g.vertT);
    Pstream::gatherList(g.vertD);
    Pstream::gatherList(g.vertC);

    return g;
}


void Foam::functionObjects::radiationDose::writeDoseCsv
(
    const AggregatedTracks& a,
    const fileName& outDir
) const
{
    OFstream os(outDir/"doseDistribution.csv");
    os  << "# trackId,endReason,time_s,dose_mJ_cm2,xEnd,yEnd,zEnd" << nl;
    forAll(a.trackId, i)
    {
        const auto er = static_cast<dose::dosePathParticle::endReason>
            (a.endReason[i]);
        const vector& x = a.xEnd[i];
        os  << a.trackId[i] << ','
            << dose::dosePathParticle::endReasonNames[er] << ','
            << a.tEnd[i] << ','
            << a.dose[i] << ','
            << x.x() << ',' << x.y() << ',' << x.z() << nl;
    }
}


void Foam::functionObjects::radiationDose::writeSummary
(
    const AggregatedTracks& a,
    const fileName& outDir
) const
{
    DynamicList<scalar> doses;
    const label totalSeeded = a.trackId.size();
    label nEscaped = 0, nTimedOut = 0, nStuck = 0, nTerm = 0, nActive = 0;
    forAll(a.trackId, i)
    {
        const auto er = static_cast<dose::dosePathParticle::endReason>
            (a.endReason[i]);
        switch (er)
        {
            case dose::dosePathParticle::endReason::escaped:
                ++nEscaped;
                doses.append(a.dose[i]);
                break;
            case dose::dosePathParticle::endReason::timedOut:
                ++nTimedOut; break;
            case dose::dosePathParticle::endReason::stuck:
                ++nStuck;    break;
            case dose::dosePathParticle::endReason::terminated:
                ++nTerm;     break;
            case dose::dosePathParticle::endReason::active:
                ++nActive;   break;
        }
    }

    const label N = doses.size();
    scalar mean = 0, minD = great, maxD = 0;
    for (label i = 0; i < N; ++i)
    {
        mean += doses[i];
        minD = min(minD, doses[i]);
        maxD = max(maxD, doses[i]);
    }
    if (N > 0) mean /= N;

    scalar var = 0;
    for (label i = 0; i < N; ++i)
    {
        var += sqr(doses[i] - mean);
    }
    if (N > 1) var /= (N - 1);
    const scalar stdev = sqrt(var);

    OFstream os(outDir/"summary.dat");
    os  << "# radiationDose summary" << nl
        << "totalSeeded     " << totalSeeded << nl
        << "escaped         " << nEscaped << nl
        << "timedOut        " << nTimedOut << nl
        << "stuck           " << nStuck << nl
        << "terminated      " << nTerm << nl
        << "stillActive     " << nActive << nl
        << "meanDose_mJcm2  " << mean << nl
        << "stdevDose_mJcm2 " << stdev << nl
        << "minDose_mJcm2   " << (N > 0 ? minD : 0) << nl
        << "maxDose_mJcm2   " << maxD << nl
        << "# Dose statistics and logReduction lines use the escaped"
        << " population only (N = " << N << ")." << nl
        << "# Tracks that timed out, got stuck, or hit maxDose are not"
        << " representative of reactor throughput and are excluded." << nl;

    forAll(kInact_, i)
    {
        const scalar k = kInact_[i];
        if (N == 0)
        {
            os  << "logReduction_k=" << k << "  0" << nl;
            continue;
        }
        scalar Ssum = 0;
        for (label j = 0; j < N; ++j)
        {
            Ssum += exp(-k*doses[j]);
        }
        const scalar S = Ssum/N;
        const scalar logR = -log10(max(S, small));
        os  << "logReduction_k=" << k << "  " << logR << nl;
    }
}


void Foam::functionObjects::radiationDose::writeVtkTrajectories
(
    const GatheredTracks& g,
    const fileName& outFile
) const
{
    // Legacy ASCII VTK PolyData. ParaView reads this directly as
    // streamlines coloured by per-vertex time / dose / cell, with
    // per-track (cell-data) trackId / endReason for filtering.
    //
    // The legacy single-file format is hand-writable without an XML
    // library and opens in ParaView identically to .vtp.
    //
    // Tracks excluded from VTK have nVerts == 0 (already filtered in
    // gatherTracks): a particle that became `stuck` before any
    // successful step has only its seed point and would form a
    // zero-length line.

    // Total counts across all ranks (already filtered: nVerts==0 for
    // skipped tracks).
    label nPoints = 0, nLines = 0, nLineConn = 0;
    forAll(g.nVerts, r)
    {
        forAll(g.nVerts[r], i)
        {
            const label nv = g.nVerts[r][i];
            if (nv == 0) continue;
            nPoints   += nv;
            nLines    += 1;
            nLineConn += nv + 1;        // [count, idx0, idx1, ...]
        }
    }

    OFstream os(outFile);
    os.precision(8);

    os  << "# vtk DataFile Version 3.0" << nl
        << "radiationDose trajectories" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << nPoints << " float" << nl;

    // Walk each rank's flattened vertex arrays in lockstep with its
    // nVerts list. The vertex offset within rank r increments by
    // nVerts[r][i] for each emitted track.
    auto walkVerts =
        [&](auto&& emit)
        {
            forAll(g.nVerts, r)
            {
                label vOff = 0;
                forAll(g.nVerts[r], i)
                {
                    const label nv = g.nVerts[r][i];
                    for (label j = 0; j < nv; ++j)
                    {
                        emit(r, vOff + j);
                    }
                    vOff += nv;
                }
            }
        };

    walkVerts
    (
        [&](label r, label v)
        {
            const vector& x = g.vertX[r][v];
            os << x.x() << ' ' << x.y() << ' ' << x.z() << nl;
        }
    );

    os  << "LINES " << nLines << ' ' << nLineConn << nl;
    {
        label idx = 0;
        forAll(g.nVerts, r)
        {
            forAll(g.nVerts[r], i)
            {
                const label nv = g.nVerts[r][i];
                if (nv == 0) continue;
                os << nv;
                for (label j = 0; j < nv; ++j)
                {
                    os << ' ' << (idx + j);
                }
                os << nl;
                idx += nv;
            }
        }
    }

    os  << "POINT_DATA " << nPoints << nl
        << "SCALARS time_s float 1" << nl
        << "LOOKUP_TABLE default" << nl;
    walkVerts
    (
        [&](label r, label v) { os << vtkSafeFloat(g.vertT[r][v]) << nl; }
    );

    os  << "SCALARS dose_mJcm2 float 1" << nl
        << "LOOKUP_TABLE default" << nl;
    walkVerts
    (
        [&](label r, label v) { os << vtkSafeFloat(g.vertD[r][v]) << nl; }
    );

    os  << "SCALARS cell int 1" << nl
        << "LOOKUP_TABLE default" << nl;
    walkVerts
    (
        [&](label r, label v) { os << g.vertC[r][v] << nl; }
    );

    os  << "CELL_DATA " << nLines << nl
        << "SCALARS trackId int 1" << nl
        << "LOOKUP_TABLE default" << nl;
    forAll(g.nVerts, r)
    {
        forAll(g.nVerts[r], i)
        {
            if (g.nVerts[r][i] == 0) continue;
            os << g.globalTrackId[r][i] << nl;
        }
    }

    // endReason as integer enum index, keyed by endReasonNames in the
    // particle header. ParaView's "Threshold" filter on this scalar
    // is the cleanest way to isolate (e.g.) only escaped tracks.
    os  << "SCALARS endReason int 1" << nl
        << "LOOKUP_TABLE default" << nl;
    forAll(g.nVerts, r)
    {
        forAll(g.nVerts[r], i)
        {
            if (g.nVerts[r][i] == 0) continue;
            os << g.endReason[r][i] << nl;
        }
    }

    // Per-track final dose. Duplicates the last per-vertex dose_mJcm2
    // value, but having it as cell data lets ParaView's "Threshold"
    // filter slice whole tracks (under/over-dosed populations) without
    // a Cell Data To Point Data conversion or a join against the CSV.
    os  << "SCALARS finalDose_mJcm2 float 1" << nl
        << "LOOKUP_TABLE default" << nl;
    forAll(g.nVerts, r)
    {
        forAll(g.nVerts[r], i)
        {
            if (g.nVerts[r][i] == 0) continue;
            os << vtkSafeFloat(g.dose[r][i]) << nl;
        }
    }
}


void Foam::functionObjects::radiationDose::writePvdWrapper
(
    const AggregatedTracks& a,
    const fileName& outDir
) const
{
    // VTK Collection (.pvd) wrapper. ParaView opens a single .pvd and
    // pulls in every referenced VTK file as one logical dataset with
    // a common Threshold target. We use one timestep with multiple
    // `part` entries (the parallel-partition semantics also serve as
    // a "load these together" hint for ParaView's reader).
    OFstream os(outDir/"trajectories.pvd");
    os  << "<?xml version=\"1.0\"?>" << nl
        << "<VTKFile type=\"Collection\" version=\"0.1\""
        << " byte_order=\"LittleEndian\">" << nl
        << "  <Collection>" << nl;
    forAll(a.vtkFile, b)
    {
        os  << "    <DataSet timestep=\"0\" part=\"" << b
            << "\" file=\"" << a.vtkFile[b] << "\"/>" << nl;
    }
    os  << "  </Collection>" << nl
        << "</VTKFile>" << nl;
}


// ************************************************************************* //
