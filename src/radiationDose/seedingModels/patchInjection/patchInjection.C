/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "patchInjection.H"
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

    // Resolve patch names to IDs and accumulate total face area
    labelList patchIDs(patchNames_.size());
    scalar totalArea = 0;
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
        const polyPatch& pp = bMesh[patchi];
        forAll(pp, fi)
        {
            totalArea += pp[fi].mag(meshPoints);
        }
    }

    if (totalArea <= 0)
    {
        FatalErrorInFunction
            << "Total area of seeding patches is zero"
            << exit(FatalError);
    }

    // Fractional particle counts per face, with stochastic rounding so that
    // expected total = nParticles_.
    DynamicList<trackPoint> seeds(nParticles_);
    label launched = 0;
    label patchIndex = 0;
    forAll(patchIDs, ip)
    {
        const polyPatch& pp = bMesh[patchIDs[ip]];
        const labelUList& cells = pp.faceCells();
        const bool last = (ip == patchIDs.size() - 1);

        forAll(pp, fi)
        {
            const face& f = pp[fi];
            const scalar A = f.mag(meshPoints);
            const scalar expected = scalar(nParticles_)*A/totalArea;
            label k = label(expected);
            const scalar frac = expected - scalar(k);
            if (rng.scalar01() < frac) { ++k; }

            // For the very last face, make up any rounding deficit so that the
            // total particle count is exactly nParticles_ on average.
            if (last && fi == pp.size() - 1)
            {
                k = max(k, nParticles_ - launched);
            }

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
                ++launched;
                ++patchIndex;
            }
        }
    }

    List<trackPoint> result;
    result.transfer(seeds);
    return result;
}


// ************************************************************************* //
