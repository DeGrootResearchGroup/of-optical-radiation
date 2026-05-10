/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "patchInjection.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(patchInjection, 0);
    addToRunTimeSelectionTable(seedingModel, patchInjection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::patchInjection::patchInjection
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    seedingModel(dict, mesh),
    patchNames_(dict.lookup("patches")),
    nParticles_(readLabel(dict.lookup("nParticles")))
{
    if (nParticles_ <= 0)
    {
        FatalErrorInFunction
            << "nParticles must be > 0, got " << nParticles_
            << exit(FatalError);
    }
    if (patchNames_.empty())
    {
        FatalErrorInFunction
            << "patchInjection requires at least one patch name"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::dose::trackPoint>
Foam::dose::patchInjection::seed(randomGenerator& rng) const
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const pointField& meshPoints = mesh_.points();

    // Resolve patch names to IDs (every rank knows the boundary
    // mesh structure even when the local patch has zero faces).
    labelList patchIDs(patchNames_.size());
    forAll(patchNames_, i)
    {
        const label patchi = bMesh.findIndex(patchNames_[i]);
        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Patch " << patchNames_[i] << " not found in mesh"
                << exit(FatalError);
        }
        patchIDs[i] = patchi;
    }

    // Local patch face area on this rank.
    scalar localArea = 0;
    forAll(patchIDs, ip)
    {
        const polyPatch& pp = bMesh[patchIDs[ip]];
        forAll(pp, fi)
        {
            localArea += pp[fi].mag(meshPoints);
        }
    }

    // Sum across ranks. In serial, reduce is a no-op and totalArea
    // == localArea. In parallel, ranks where the patch has zero
    // local faces contribute 0 and proceed without seeding any
    // particles below.
    scalar totalArea = localArea;
    reduce(totalArea, sumOp<scalar>());

    if (totalArea <= 0)
    {
        FatalErrorInFunction
            << "Total area of seeding patches is zero across all ranks"
            << exit(FatalError);
    }

    // Per-face fractional particle count weighted by the GLOBAL area,
    // with stochastic rounding so the expected total across all
    // ranks is nParticles_. Variance is O(sqrt(nParticles_)). The
    // pre-parallel implementation also added a deterministic
    // "deficit correction" on the last face to land exactly on
    // nParticles_ in serial; that depended on each rank knowing the
    // running global launched count, which is not free in parallel.
    // We accept the small extra variance in exchange for a clean
    // global-area-driven scheme that works the same way serial and
    // parallel.
    DynamicList<trackPoint> seeds;
    forAll(patchIDs, ip)
    {
        const polyPatch& pp = bMesh[patchIDs[ip]];
        const labelUList& cells = pp.faceCells();

        forAll(pp, fi)
        {
            const face& f = pp[fi];
            const scalar A = f.mag(meshPoints);
            const scalar expected = scalar(nParticles_)*A/totalArea;
            label k = label(expected);
            const scalar frac = expected - scalar(k);
            if (rng.scalar01() < frac) { ++k; }
            if (k <= 0) continue;

            // Triangulate the face fan-style around its centroid.
            const point fc = f.centre(meshPoints);
            const label nv = f.size();

            // Sub-triangle areas for weighted sampling
            scalarList triA(nv);
            scalar sumA = 0;
            for (label v = 0; v < nv; ++v)
            {
                const point& a = meshPoints[f[v]];
                const point& b = meshPoints[f[(v + 1) % nv]];
                triA[v] = 0.5*mag((a - fc) ^ (b - fc));
                sumA += triA[v];
            }

            for (label n = 0; n < k; ++n)
            {
                // Pick a sub-triangle weighted by area
                scalar r = rng.scalar01()*sumA;
                label tri = 0;
                scalar acc = triA[0];
                while (r > acc && tri < nv - 1)
                {
                    ++tri;
                    acc += triA[tri];
                }

                // Sample uniformly within the chosen triangle
                scalar u = rng.scalar01();
                scalar v = rng.scalar01();
                if (u + v > 1) { u = 1 - u; v = 1 - v; }

                const point& a = meshPoints[f[tri]];
                const point& b = meshPoints[f[(tri + 1) % nv]];
                const point p = fc + u*(a - fc) + v*(b - fc);

                seeds.append(trackPoint(p, 0, 0, cells[fi]));
            }
        }
    }

    List<trackPoint> result;
    result.transfer(seeds);
    return result;
}


// ************************************************************************* //
