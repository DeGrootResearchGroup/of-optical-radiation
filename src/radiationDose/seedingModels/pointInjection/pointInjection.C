/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "pointInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(pointInjection, 0);
    addToRunTimeSelectionTable(seedingModel, pointInjection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::pointInjection::pointInjection
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    seedingModel(dict, mesh),
    nParticles_(readLabel(dict.lookup("nParticles"))),
    shape_(shapeType::sphere),
    centre_(Zero),
    radius_(0),
    size_(Zero),
    maxAttempts_(dict.lookupOrDefault<label>("maxAttempts", 1000))
{
    if (nParticles_ <= 0)
    {
        FatalErrorInFunction
            << "nParticles must be > 0, got " << nParticles_
            << exit(FatalError);
    }

    const dictionary& regionDict = dict.subDict("region");
    centre_ = point(regionDict.lookup("centre"));

    const word shapeName(regionDict.lookup("type"));
    if (shapeName == "sphere")
    {
        shape_ = shapeType::sphere;
        radius_ = readScalar(regionDict.lookup("radius"));
        if (radius_ <= 0)
        {
            FatalErrorInFunction
                << "sphere radius must be > 0, got " << radius_
                << exit(FatalError);
        }
    }
    else if (shapeName == "box")
    {
        shape_ = shapeType::box;
        size_ = vector(regionDict.lookup("size"));
        if (size_.x() <= 0 || size_.y() <= 0 || size_.z() <= 0)
        {
            FatalErrorInFunction
                << "box size components must all be > 0, got " << size_
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown region type " << shapeName
            << "; expected 'sphere' or 'box'"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::dose::trackPoint>
Foam::dose::pointInjection::seed(randomGenerator& rng) const
{
    DynamicList<trackPoint> seeds(nParticles_);

    // The mesh-cached meshSearch handles the octree-backed point
    // location used for rejecting candidates that lie outside this
    // rank's local mesh portion (and outside the global mesh in
    // serial). One construction; subsequent New() calls are O(1).
    const meshSearch& ms = meshSearch::New(mesh_);

    // Half-extents of the rejection bounding box. For a sphere this is
    // the radius along all three axes; for a box it's size/2.
    const vector halfExtent =
        (shape_ == shapeType::sphere)
      ? vector(radius_, radius_, radius_)
      : 0.5*size_;

    // The inner draw loop keeps every rank's RNG in lockstep: the only
    // randomness is the candidate-point generation and the shape test;
    // both are deterministic functions of the RNG state, so every rank
    // arrives at the same accepted-shape candidate after the same number
    // of draws. Only mesh_.findCell() differs per rank (each rank's
    // mesh is a different subdomain), so exactly one rank accepts each
    // particle and the global total stays bounded by nParticles_.
    for (label n = 0; n < nParticles_; ++n)
    {
        bool haveCandidate = false;
        point x(Zero);

        for (label attempt = 0; attempt < maxAttempts_ && !haveCandidate; ++attempt)
        {
            // Uniform in [-1, 1]^3 — three RNG draws per attempt
            const vector u
            (
                2*rng.scalar01() - 1,
                2*rng.scalar01() - 1,
                2*rng.scalar01() - 1
            );
            x = centre_ + cmptMultiply(halfExtent, u);

            if (shape_ == shapeType::sphere)
            {
                // |u|^2 < 1 iff x is inside the sphere of radius
                // halfExtent; with halfExtent = (R, R, R) this is
                // exactly |x - centre| < R.
                if (magSqr(u) > 1) continue;
            }
            // box: every cube draw is already inside the box.
            haveCandidate = true;
        }

        if (!haveCandidate)
        {
            // maxAttempts_ exceeded purely by shape rejection - implies
            // the rejection rate is far worse than ~50 %, which only
            // happens if the user's region degenerated. Bail out.
            FatalErrorInFunction
                << "pointInjection: could not draw a shape-accepted "
                << "candidate in " << maxAttempts_ << " attempts. "
                << "Check region configuration."
                << exit(FatalError);
        }

        // Rank-specific: each rank only owns part of the mesh. A point
        // in another rank's subdomain returns -1 here and is silently
        // skipped on this rank; the matching rank will accept it.
        // A point outside the global domain is rejected on every rank
        // (-> particle count comes in below nParticles_).
        const label celli = ms.findCell(x);
        if (celli >= 0)
        {
            seeds.append(trackPoint(x, 0, 0, celli));
        }
    }

    List<trackPoint> result;
    result.transfer(seeds);
    return result;
}


// ************************************************************************* //
