/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "radiationDose.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"

#ifdef _OPENMP
    #include <omp.h>
#else
    static inline int omp_get_thread_num() { return 0; }
    static inline int omp_get_num_threads() { return 1; }
#endif

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
    maxStepsPerTrack_(100000)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

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
    dispersion_ = dose::dispersionModel::New(dict.subDict("dispersion"), mesh_);

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
    maxStepsPerTrack_ =
        intDict.lookupOrDefault<label>("maxStepsPerTrack", 100000);

    const dictionary& outDict = dict.subDict("output");
    kInact_ = outDict.lookupOrDefault<scalarList>("kInact", scalarList());

    return true;
}


Foam::wordList Foam::functionObjects::radiationDose::fields() const
{
    DynamicList<word> f;
    f.append(UName_);
    f.append(GName_);
    if (dispersion_.valid())
    {
        const wordList extra = dispersion_->requiredFields();
        forAll(extra, i) { f.append(extra[i]); }
    }
    return wordList(f);
}


bool Foam::functionObjects::radiationDose::execute()
{
    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& G =
        mesh_.lookupObject<volScalarField>(GName_);

    randomGenerator rng{randomGenerator::seed(randomSeed_)};

    autoPtr<interpolation<vector>> UintPtr =
        interpolation<vector>::New("cellPoint", U);
    autoPtr<interpolation<scalar>> GintPtr =
        interpolation<scalar>::New("cellPoint", G);

    Info<< type() << ": seeding particles..." << endl;
    List<dose::trackPoint> seeds = seeding_->seed(rng);
    Info<< type() << ": seeded " << seeds.size() << " particles" << endl;

    tracks_.clear();
    tracks_.setSize(seeds.size());
    forAll(seeds, i)
    {
        tracks_.set(i, new dose::track(i));
        tracks_[i].append(seeds[i]);
        // Per-track dispersion state lives on the track, not on the
        // model — this is what makes the integration loop thread-safe.
        tracks_[i].setDispState(dispersion_->newState());
    }

    // Threaded integration. Each particle's trajectory is independent
    // of every other particle's (mesh + fields are read-only, dispersion
    // state is per-track), so we partition tracks across threads and
    // give each thread its own RNG (seeded deterministically from the
    // master seed + thread id, so runs are reproducible).
    #ifdef _OPENMP
    Info<< type() << ": integrating with " << omp_get_max_threads()
        << " thread(s)..." << endl;
    #else
    Info<< type() << ": integrating (serial)..." << endl;
    #endif

    #pragma omp parallel
    {
        randomGenerator threadRng
        {
            randomGenerator::seed(randomSeed_ + omp_get_thread_num())
        };

        #pragma omp for schedule(dynamic, 50)
        for (label i = 0; i < tracks_.size(); ++i)
        {
            integrateTrack(tracks_[i], UintPtr(), GintPtr(), threadRng);
        }
    }

    label nEscaped = 0, nTimedOut = 0, nStuck = 0, nLeft = 0, nOther = 0;
    forAll(tracks_, i)
    {
        switch (tracks_[i].end())
        {
            case dose::track::endReason::escaped:    ++nEscaped;  break;
            case dose::track::endReason::timedOut:   ++nTimedOut; break;
            case dose::track::endReason::stuck:      ++nStuck;    break;
            case dose::track::endReason::leftDomain: ++nLeft;     break;
            default:                                 ++nOther;    break;
        }
    }
    Info<< type() << ": done. "
        << nEscaped << " escaped, "
        << nTimedOut << " timed out, "
        << nStuck << " stuck, "
        << nLeft << " left domain, "
        << nOther << " other" << endl;

    return true;
}


void Foam::functionObjects::radiationDose::integrateTrack
(
    dose::track& tr,
    const interpolation<vector>& Uint,
    const interpolation<scalar>& Gint,
    randomGenerator& rng
) const
{
    // trackToFace integration:
    //
    //   Each outer iteration samples a velocity V at the current position
    //   and walks face-by-face through the dt budget. The inner loop
    //   advances cell-by-cell across internal faces (updating celli via
    //   faceOwner/faceNeighbour) and, on hitting a boundary face, either
    //   terminates (escape / stuck) or specularly reflects V and continues
    //   inside the same cell. Cell index is maintained explicitly through
    //   face crossings so no findCell calls are needed mid-track.
    //
    //   This eliminates the v0.1 leftDomain class of bug, which arose from
    //   the post-step findCell + boundary-of-celli_prev scan missing
    //   segment-face intersections on curved walls.
    //
    //   Notation: for each face fi of celli the area-weighted normal Sf[fi]
    //   points OUT of the face's owner cell. Whether celli is the owner or
    //   the neighbour determines the "outward" sign convention for that face.

    const scalarField& cellVolumes = mesh_.V();
    const cellList& meshCells = mesh_.cells();
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const label nInternal = mesh_.nInternalFaces();
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    // Cap reflections within a single outer step: catches pathological
    // corner cases where a particle would bounce indefinitely. Set
    // generously: at high turbulent intensity, DRW fluctuations can push
    // a near-wall particle back into the wall repeatedly across many
    // velocity samples; we want to allow the natural physics rather than
    // class genuinely-tracked particles as stuck. The cap exists only
    // to bound runtime in pathological corners.
    constexpr label maxReflectionsPerStep = 100;

    label step = 0;
    while (tr.end() == dose::track::endReason::active
        && step < maxStepsPerTrack_)
    {
        const dose::trackPoint p = tr.back();
        if (p.celli < 0)
        {
            tr.setEnd(dose::track::endReason::leftDomain);
            break;
        }

        const vector Umean = Uint.interpolate(p.x, p.celli);
        vector V = Umean
            + dispersion_->fluctuation
              (
                  tr.dispState(),
                  p.x,
                  p.celli,
                  dtMax_,
                  rng
              );

        if (mag(V) < small)
        {
            tr.setEnd(dose::track::endReason::stuck);
            break;
        }

        // CFL-bounded budget for this outer iteration
        const scalar cellSize = cbrt(cellVolumes[p.celli]);
        scalar dtRemaining = min(dtMax_, cflMax_*cellSize/mag(V));

        point x = p.x;
        label celli = p.celli;
        scalar D = p.D;
        scalar t_now = p.t;
        // Clamp interpolated G to [0, +inf): cellPoint interpolation can
        // produce small negative values near boundaries on cells where the
        // tet-decomposition's vertex weights overshoot, even when the
        // underlying field is strictly non-negative. Clamping here keeps the
        // dose accumulation monotonic.
        scalar G_here = max(scalar(0), Gint.interpolate(x, celli));
        label nReflected = 0;

        bool terminated = false;
        while (dtRemaining > small)
        {
            // Find first face crossing
            scalar tHit = dtRemaining;
            label hitFi = -1;

            const cell& c = meshCells[celli];
            forAll(c, i)
            {
                const label fi = c[i];

                // outward = Sf when celli is owner, -Sf when celli is neighbour
                const scalar sign = (faceOwner[fi] == celli) ? 1.0 : -1.0;
                const scalar denom = sign*(Sf[fi] & V);
                if (denom <= small)
                {
                    // V doesn't point out through this face
                    continue;
                }

                const scalar t = sign*((Cf[fi] - x) & Sf[fi])/denom;
                if (t < -small)
                {
                    continue;
                }

                if (t < tHit)
                {
                    tHit = t;
                    hitFi = fi;
                }
            }

            // Advance to the hit (or to dtRemaining if none) and accrue dose
            const scalar tStep = max(scalar(0), tHit);
            const point xNext = x + tStep*V;
            const scalar G_at_next =
                max(scalar(0), Gint.interpolate(xNext, celli));
            D += 0.5*(G_here + G_at_next)*tStep*Wm2_s_to_mJcm2;
            t_now += tStep;
            x = xNext;
            G_here = G_at_next;
            dtRemaining -= tStep;

            if (hitFi < 0)
            {
                // No face reachable in remaining dt; particle stays in cell
                break;
            }

            if (hitFi < nInternal)
            {
                // Internal face: walk to neighbour
                celli = (faceOwner[hitFi] == celli)
                    ? faceNeighbour[hitFi] : faceOwner[hitFi];
                G_here = max(scalar(0), Gint.interpolate(x, celli));
                continue;
            }

            // Boundary face
            const label patchi = bMesh.whichPatch(hitFi);
            if (escapePatchIDs_.found(patchi))
            {
                tr.append(dose::trackPoint(x, t_now, D, celli));
                tr.setEnd(dose::track::endReason::escaped);
                terminated = true;
                break;
            }
            if (!wallReflection_)
            {
                tr.append(dose::trackPoint(x, t_now, D, celli));
                tr.setEnd(dose::track::endReason::stuck);
                terminated = true;
                break;
            }

            // Reflect V about the face normal and stay in the same cell.
            const scalar Smag = mag(Sf[hitFi]);
            if (Smag < small)
            {
                tr.append(dose::trackPoint(x, t_now, D, celli));
                tr.setEnd(dose::track::endReason::stuck);
                terminated = true;
                break;
            }
            const vector n = Sf[hitFi]/Smag;
            V = V - 2.0*(V & n)*n;

            ++nReflected;
            if (nReflected > maxReflectionsPerStep)
            {
                tr.append(dose::trackPoint(x, t_now, D, celli));
                tr.setEnd(dose::track::endReason::stuck);
                terminated = true;
                break;
            }
        }

        if (terminated)
        {
            break;
        }

        tr.append(dose::trackPoint(x, t_now, D, celli));

        if (maxTime_ > 0 && t_now >= maxTime_)
        {
            tr.setEnd(dose::track::endReason::timedOut);
            break;
        }
        if (maxDose_ > 0 && D >= maxDose_)
        {
            tr.setEnd(dose::track::endReason::terminated);
            break;
        }

        ++step;
    }

    if (step >= maxStepsPerTrack_
     && tr.end() == dose::track::endReason::active)
    {
        tr.setEnd(dose::track::endReason::stuck);
    }
}


bool Foam::functionObjects::radiationDose::write()
{
    writeDoseCsv();
    writeSummary();
    return true;
}


void Foam::functionObjects::radiationDose::writeDoseCsv() const
{
    const fileName outDir =
        time_.path()/"postProcessing"/name()/time_.name();
    mkDir(outDir);

    OFstream os(outDir/"doseDistribution.csv");
    os  << "# trackId,endReason,time_s,dose_mJ_cm2,xEnd,yEnd,zEnd" << nl;
    forAll(tracks_, i)
    {
        const dose::track& tr = tracks_[i];
        const dose::trackPoint& p = tr.back();
        os  << tr.id() << ','
            << dose::track::endReasonNames[tr.end()] << ','
            << p.t << ','
            << p.D << ','
            << p.x.x() << ',' << p.x.y() << ',' << p.x.z() << nl;
    }
}


void Foam::functionObjects::radiationDose::writeSummary() const
{
    const fileName outDir =
        time_.path()/"postProcessing"/name()/time_.name();
    mkDir(outDir);

    // Count active (escaped) tracks and gather their final doses
    DynamicList<scalar> doses(tracks_.size());
    forAll(tracks_, i)
    {
        if (tracks_[i].end() == dose::track::endReason::escaped)
        {
            doses.append(tracks_[i].finalDose());
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
        << "totalSeeded     " << tracks_.size() << nl
        << "escaped         " << N << nl
        << "meanDose_mJcm2  " << mean << nl
        << "stdevDose_mJcm2 " << stdev << nl
        << "minDose_mJcm2   " << (N > 0 ? minD : 0) << nl
        << "maxDose_mJcm2   " << maxD << nl;

    // Log reduction = -log10( (1/N) sum_i exp(-k * D_i) )
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


// ************************************************************************* //
