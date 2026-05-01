/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "radiationDose.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"
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
    return wordList({UName_, GName_});
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

    const meshSearch& search = meshSearch::New(mesh_);

    Info<< type() << ": seeding particles..." << endl;
    List<dose::trackPoint> seeds = seeding_->seed(rng);
    Info<< type() << ": seeded " << seeds.size() << " particles" << endl;

    tracks_.clear();
    tracks_.setSize(seeds.size());
    forAll(seeds, i)
    {
        tracks_.set(i, new dose::track(i));
        tracks_[i].append(seeds[i]);
    }

    Info<< type() << ": integrating..." << endl;
    forAll(tracks_, i)
    {
        integrateTrack(tracks_[i], UintPtr(), GintPtr(), search, rng);
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
    const meshSearch& search,
    randomGenerator& rng
) const
{
    const scalarField& cellVolumes = mesh_.V();

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
        const vector Ufluc =
            dispersion_->fluctuation(p.x, p.celli, dtMax_, rng);
        const vector V = Umean + Ufluc;

        const scalar Vmag = mag(V);
        if (Vmag < small)
        {
            tr.setEnd(dose::track::endReason::stuck);
            break;
        }

        const scalar cellSize = cbrt(cellVolumes[p.celli]);
        scalar dt = min(dtMax_, cflMax_*cellSize/Vmag);

        // RK2 midpoint
        point x_mid = p.x + 0.5*dt*V;
        label celli_mid = search.findCell(x_mid);
        if (celli_mid < 0)
        {
            // Mid-step crossed boundary; halve dt and retry once
            dt *= 0.5;
            x_mid = p.x + 0.5*dt*V;
            celli_mid = search.findCell(x_mid);
        }

        vector V_mid;
        if (celli_mid >= 0)
        {
            const vector Umean_mid = Uint.interpolate(x_mid, celli_mid);
            const vector Ufluc_mid =
                dispersion_->fluctuation(x_mid, celli_mid, 0.5*dt, rng);
            V_mid = Umean_mid + Ufluc_mid;
        }
        else
        {
            // Fallback to forward Euler if mid-point still outside
            V_mid = V;
        }

        const point x_new = p.x + dt*V_mid;
        const label celli_new = search.findCell(x_new);

        if (celli_new < 0)
        {
            handleBoundaryHit(tr, p.x, p.celli, V_mid, dt, Gint);
            break;
        }

        // Accumulate dose using trapezoidal rule on G
        const scalar G_n = Gint.interpolate(p.x, p.celli);
        const scalar G_new = Gint.interpolate(x_new, celli_new);
        const scalar G_avg = 0.5*(G_n + G_new);
        const scalar D_new = p.D + G_avg*dt*Wm2_s_to_mJcm2;
        const scalar t_new = p.t + dt;

        tr.append(dose::trackPoint(x_new, t_new, D_new, celli_new));

        if (maxTime_ > 0 && t_new >= maxTime_)
        {
            tr.setEnd(dose::track::endReason::timedOut);
            break;
        }
        if (maxDose_ > 0 && D_new >= maxDose_)
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


void Foam::functionObjects::radiationDose::handleBoundaryHit
(
    dose::track& tr,
    const point& x_prev,
    label celli_prev,
    const vector& V,
    scalar dt,
    const interpolation<scalar>& Gint
) const
{
    // Find which boundary face of celli_prev the segment x_prev -> x_prev+dt*V
    // crosses first.
    const cell& c = mesh_.cells()[celli_prev];
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Cf = mesh_.faceCentres();
    const label nInternal = mesh_.nInternalFaces();

    scalar tHit = great;
    label hitFacei = -1;
    forAll(c, i)
    {
        const label fi = c[i];
        if (fi < nInternal) continue;

        const scalar denom = Sf[fi] & V;
        if (denom <= 0) continue;  // not exiting through this face

        const scalar tCross = ((Cf[fi] - x_prev) & Sf[fi]) / (denom*dt);
        if (tCross < 0 || tCross > 1.0 + small) continue;

        if (tCross < tHit)
        {
            tHit = tCross;
            hitFacei = fi;
        }
    }

    dose::track::endReason reason = dose::track::endReason::leftDomain;
    point xHit = x_prev + dt*V;     // fallback if no face found
    scalar dtPartial = dt;

    if (hitFacei >= 0)
    {
        dtPartial = max(small, tHit*dt);
        xHit = x_prev + dtPartial*V;
        const label hitPatchi =
            mesh_.boundaryMesh().whichPatch(hitFacei);
        reason = escapePatchIDs_.found(hitPatchi)
            ? dose::track::endReason::escaped
            : dose::track::endReason::stuck;
    }

    // Trapezoidal dose contribution over the partial step (G evaluated in
    // the previous cell at both endpoints, since the hit point lies on
    // celli_prev's boundary)
    const scalar G_n = Gint.interpolate(x_prev, celli_prev);
    const scalar G_hit = Gint.interpolate(xHit, celli_prev);
    const scalar G_avg = 0.5*(G_n + G_hit);
    const scalar D_new = tr.back().D + G_avg*dtPartial*Wm2_s_to_mJcm2;
    const scalar t_new = tr.back().t + dtPartial;

    tr.append(dose::trackPoint(xHit, t_new, D_new, celli_prev));
    tr.setEnd(reason);
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
