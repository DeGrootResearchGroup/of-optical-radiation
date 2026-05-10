/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "radiationDose.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolationCellPoint.H"
#include "meshSearch.H"
#include "polyBoundaryMesh.H"

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
    UName_("U"),
    GName_("G"),
    randomSeed_(0),
    maxTime_(0),
    maxDose_(0),
    wallReflection_(true),
    dtMax_(0.01),
    cflMax_(0.5),
    maxOuterSteps_(100000),
    writeVtk_(true)
{
    read(dict);
}


Foam::functionObjects::radiationDose::~radiationDose()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::radiationDose::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

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
        const label patchi =
            mesh_.boundaryMesh().findIndex(escapePatches[i]);
        if (patchi < 0)
        {
            FatalErrorInFunction
                << "escapePatch " << escapePatches[i]
                << " not found in mesh"
                << exit(FatalError);
        }
        escapePatchIDs_.insert(patchi);
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
    writeVtk_ = outDict.lookupOrDefault<Switch>("writeVtk", true);

    return true;
}


Foam::wordList Foam::functionObjects::radiationDose::fields() const
{
    DynamicList<word> f;
    f.append(UName_);
    f.append(GName_);
    // Probe the dispersion and motion models for any extra fields they
    // need (e.g. k, epsilon for DRW). Done with throw-away instances —
    // the models themselves are rebuilt fresh in execute().
    autoPtr<dose::dispersionModel> dispProbe =
        dose::dispersionModel::New(dispersionDict_, mesh_);
    const wordList dispExtra = dispProbe->requiredFields();
    forAll(dispExtra, i) { f.append(dispExtra[i]); }

    autoPtr<dose::motionModel> motProbe =
        dose::motionModel::New(motionDict_, mesh_);
    const wordList motExtra = motProbe->requiredFields();
    forAll(motExtra, i) { f.append(motExtra[i]); }

    return wordList(f);
}


bool Foam::functionObjects::radiationDose::execute()
{
    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& G =
        mesh_.lookupObject<volScalarField>(GName_);

    randomGenerator rng{randomGenerator::seed(randomSeed_)};

    // Build the dispersion and motion models fresh per execute() so
    // any internal counters / RNG state are reproducible.
    autoPtr<dose::dispersionModel> dispersion =
        dose::dispersionModel::New(dispersionDict_, mesh_);
    autoPtr<dose::motionModel> motion =
        dose::motionModel::New(motionDict_, mesh_);

    // Construct the cloud (or rebuild if execute() is called twice in
    // one run). The cloud takes ownership of both models.
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
            std::move(dispersion),
            std::move(motion)
        )
    );

    // Seed the cloud. The seedingModel returns initial trackPoints
    // (position + cell index); we use those to construct dosePath-
    // particles via OpenFOAM's locate-by-tree constructor.
    Info<< type() << ": seeding particles..." << endl;
    List<dose::trackPoint> seeds = seeding_->seed(rng);
    Info<< type() << ": seeded " << seeds.size() << " particles" << endl;

    const meshSearch& ms = meshSearch::New(mesh_);
    label nLocateBoundaryHits = 0;
    forAll(seeds, i)
    {
        dose::dosePathParticle* p =
            new dose::dosePathParticle
            (
                ms,
                seeds[i].x,
                seeds[i].celli,
                nLocateBoundaryHits
            );
        // Per-track dispersion and motion states are owned by the
        // particle; the cloud's models are the factories.
        p->setDispState(cloud_->dispersion().newState());
        p->setMotionState(cloud_->motion().newState());
        // Initial trajectory vertex (so the writer reports the seed
        // position even if the particle never advances). Skipped
        // when storeTrack is off — see dosePathParticle::move() for
        // the per-step variant of this same guard.
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
        cloud_->runToCompletion
        (
            UInterp,
            GInterp,
            rng,
            maxOuterSteps_
        );

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

    return true;
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
};


bool Foam::functionObjects::radiationDose::write()
{
    if (!cloud_.valid())
    {
        return true;
    }

    const GatheredTracks g = gatherTracks();

    if (Pstream::master())
    {
        // globalPath() strips the processorN suffix in parallel runs;
        // in serial it's identical to path(). The summary / CSV / VTK
        // are already global (rank 0 holds all ranks' data after the
        // gather), so they belong in the case-root postProcessing
        // tree, not under any processor directory.
        const fileName outDir =
            time_.globalPath()/"postProcessing"/name()/time_.name();
        mkDir(outDir);

        writeDoseCsv(g, outDir);
        writeSummary(g, outDir);
        if (writeVtk_)
        {
            writeVtkTrajectories(g, outDir);
        }
    }
    return true;
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
        const label nv =
            (writeVtk_ && pts.size() >= 2) ? pts.size() : 0;
        g.nVerts[me][k] = nv;
        for (label j = 0; j < nv; ++j)
        {
            vx.append(pts[j].x);
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
    const GatheredTracks& g,
    const fileName& outDir
) const
{
    OFstream os(outDir/"doseDistribution.csv");
    os  << "# trackId,endReason,time_s,dose_mJ_cm2,xEnd,yEnd,zEnd" << nl;
    forAll(g.origId, r)
    {
        forAll(g.origId[r], i)
        {
            const auto er = static_cast<dose::dosePathParticle::endReason>
                (g.endReason[r][i]);
            const vector& x = g.xEnd[r][i];
            os  << g.origProc[r][i] << '.' << g.origId[r][i] << ','
                << dose::dosePathParticle::endReasonNames[er] << ','
                << g.tEnd[r][i] << ','
                << g.dose[r][i] << ','
                << x.x() << ',' << x.y() << ',' << x.z() << nl;
        }
    }
}


void Foam::functionObjects::radiationDose::writeSummary
(
    const GatheredTracks& g,
    const fileName& outDir
) const
{
    DynamicList<scalar> doses;
    label totalSeeded = 0;
    forAll(g.origId, r)
    {
        totalSeeded += g.origId[r].size();
        forAll(g.origId[r], i)
        {
            if
            (
                static_cast<dose::dosePathParticle::endReason>
                    (g.endReason[r][i])
             == dose::dosePathParticle::endReason::escaped
            )
            {
                doses.append(g.dose[r][i]);
            }
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
        << "escaped         " << N << nl
        << "meanDose_mJcm2  " << mean << nl
        << "stdevDose_mJcm2 " << stdev << nl
        << "minDose_mJcm2   " << (N > 0 ? minD : 0) << nl
        << "maxDose_mJcm2   " << maxD << nl;

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
    const fileName& outDir
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

    OFstream os(outDir/"trajectories.vtk");
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
        [&](label r, label v) { os << g.vertT[r][v] << nl; }
    );

    os  << "SCALARS dose_mJcm2 float 1" << nl
        << "LOOKUP_TABLE default" << nl;
    walkVerts
    (
        [&](label r, label v) { os << g.vertD[r][v] << nl; }
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
            os << g.origId[r][i] << nl;
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
}


// ************************************************************************* //
