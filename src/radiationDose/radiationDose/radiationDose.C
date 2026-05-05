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
    maxOuterSteps_(100000)
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

    return true;
}


Foam::wordList Foam::functionObjects::radiationDose::fields() const
{
    DynamicList<word> f;
    f.append(UName_);
    f.append(GName_);
    // The dispersion model may need additional fields (e.g. k, epsilon
    // for DRW). Probe its requiredFields() with a temporary instance.
    autoPtr<dose::dispersionModel> probe =
        dose::dispersionModel::New(dispersionDict_, mesh_);
    const wordList extra = probe->requiredFields();
    forAll(extra, i) { f.append(extra[i]); }
    return wordList(f);
}


bool Foam::functionObjects::radiationDose::execute()
{
    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);
    const volScalarField& G =
        mesh_.lookupObject<volScalarField>(GName_);

    randomGenerator rng{randomGenerator::seed(randomSeed_)};

    // Build the dispersion model fresh per execute() so its internal
    // counter / RNG state is reproducible.
    autoPtr<dose::dispersionModel> dispersion =
        dose::dispersionModel::New(dispersionDict_, mesh_);

    // Construct the cloud (or rebuild if execute() is called twice in
    // one run). The cloud takes ownership of the dispersion model.
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
            std::move(dispersion)
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
        // Per-track dispersion state is owned by the particle; the
        // cloud's dispersion model is the factory.
        p->setDispState(cloud_->dispersion().newState());
        // Initial trajectory vertex (so the writer reports the seed
        // position even if the particle never advances).
        p->appendPoint(seeds[i]);
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


bool Foam::functionObjects::radiationDose::write()
{
    if (!cloud_.valid())
    {
        return true;
    }
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
    forAllConstIter
    (
        typename lagrangian::Cloud<dose::dosePathParticle>,
        *cloud_,
        iter
    )
    {
        const dose::dosePathParticle& p = iter();
        const vector x = p.position(mesh_);
        os  << p.origId() << ','
            << dose::dosePathParticle::endReasonNames[p.end()] << ','
            << p.t() << ','
            << p.D() << ','
            << x.x() << ',' << x.y() << ',' << x.z() << nl;
    }
}


void Foam::functionObjects::radiationDose::writeSummary() const
{
    const fileName outDir =
        time_.path()/"postProcessing"/name()/time_.name();
    mkDir(outDir);

    DynamicList<scalar> doses(cloud_->size());
    forAllConstIter
    (
        typename lagrangian::Cloud<dose::dosePathParticle>,
        *cloud_,
        iter
    )
    {
        if (iter().end() == dose::dosePathParticle::endReason::escaped)
        {
            doses.append(iter().D());
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
        << "totalSeeded     " << cloud_->size() << nl
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


// ************************************************************************* //
