/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::dosePathCloud

Description
    Cloud of dose-tracking particles, derived from Foam's
    lagrangian::Cloud<dosePathParticle>. Holds the case configuration
    that the per-particle move() needs to consult: outer-step duration
    dtMax, CFL bound cflMax, escape-patch ID set, soft-termination
    bounds maxTime / maxDose, wall-reflection switch, and the
    dispersion model.

    The standard OpenFOAM Cloud::move() advances every particle by one
    OF time step (which here is interpreted as one outer step of
    dtMax). To run particles to completion in a postProcess context
    we wrap that in runToCompletion(), which drives Cloud::move() in a
    loop until no particles are still active or a hard step cap is
    reached.

\*---------------------------------------------------------------------------*/

#ifndef dose_dosePathCloud_H
#define dose_dosePathCloud_H

#include "Cloud.H"
#include "HashSet.H"
#include "Switch.H"
#include "autoPtr.H"
#include "dispersionModel.H"
#include "dosePathParticle.H"
#include "fvMesh.H"
#include "interpolationCellPoint.H"
#include "randomGenerator.H"

#include <vector>

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class dosePathCloud Declaration
\*---------------------------------------------------------------------------*/

class dosePathCloud
:
    public lagrangian::Cloud<dosePathParticle>
{
    // Private Data — case configuration

        const fvMesh& mesh_;

        scalar dtMax_;
        scalar cflMax_;
        labelHashSet escapePatchIDs_;
        scalar maxTime_;
        scalar maxDose_;
        Switch wallReflection_;

        //- Whether dosePathParticle::move() should append a vertex to
        //  its `points_` list at the end of every outer step. Off
        //  for runs that disable VTK output (the only consumer of
        //  the stored trajectory) — bounds memory at O(N_particles)
        //  instead of O(N_particles × residence_time / dtMax).
        Switch storeTrack_;

        //- Dispersion model (RTS-selected). Lives here so the per-
        //  particle move() can sample fluctuations without needing a
        //  separate handle.
        autoPtr<dispersionModel> dispersion_;


public:

    //- Runtime type information
    TypeName("dosePathCloud");


    // Constructors

        //- Construct from a mesh and case parameters. The dispersion
        //  model is owned here (passed by autoPtr move) so its
        //  newState() factory can hand each seeded particle its own
        //  per-track state.
        dosePathCloud
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
            autoPtr<dispersionModel> dispersion
        );

        //- Disallow default bitwise copy construction
        dosePathCloud(const dosePathCloud&) = delete;


    // Member Functions

        // Access

            inline const fvMesh& mesh() const { return mesh_; }
            inline scalar dtMax() const { return dtMax_; }
            inline scalar cflMax() const { return cflMax_; }
            inline const labelHashSet& escapePatchIDs() const
            {
                return escapePatchIDs_;
            }
            inline scalar maxTime() const { return maxTime_; }
            inline scalar maxDose() const { return maxDose_; }
            inline Switch wallReflection() const { return wallReflection_; }
            inline Switch storeTrack() const { return storeTrack_; }
            inline dispersionModel& dispersion() { return *dispersion_; }
            inline const dispersionModel& dispersion() const
            {
                return *dispersion_;
            }


        // Edit

            //- Drive an outer-step loop until every particle is no
            //  longer active (escaped / timedOut / stuck /
            //  terminated) or maxOuterSteps is reached as a safety
            //  net. Returns the number of outer steps actually run.
            //
            //  In single-rank runs the inner per-particle iteration
            //  is OpenMP-parallelised across the local cloud (see
            //  moveOmp); in multi-rank MPI runs we fall back to OF's
            //  serial Cloud::move(), which carries the parallel
            //  particle handoff machinery needed at processor patches.
            //  OMP-within-MPI is a future optimisation — the per-rank
            //  transfer queues Cloud::move() builds at processor
            //  patches would need explicit thread-safe synchronisation
            //  to be safe under threading, and there is no real-world
            //  driver for it yet.
            label runToCompletion
            (
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& GInterp,
                randomGenerator& parentRng,
                const label maxOuterSteps
            );

            //- Number of particles still in the active state
            //  (endReason == active). Counted via a forAllConstIter.
            label nActive() const;


    private:

        // Private Member Functions

            //- Single outer step, OMP-parallelised across local
            //  particles. Takes pre-built per-thread trackingData
            //  and the particle pointer snapshot — `runToCompletion`
            //  hoists this state out of the outer-step loop so it
            //  is not rebuilt on every call. Each thread is bound
            //  to its own randomGenerator (seeded deterministically
            //  from the parent RNG and the thread index) via the
            //  trackingData reference; field interpolators are
            //  shared read-only.
            //
            //  Reproducibility: results are bit-for-bit reproducible
            //  for a fixed (seed, nThreads, particle ordering). The
            //  static schedule + parent-derived per-thread seeding
            //  mean particle `i` is always processed by thread
            //  (i / chunk) with the same RNG sequence — but
            //  switching nThreads reshuffles which particles see
            //  which RNG samples and therefore the per-particle
            //  dose outcomes (the ensemble statistics still converge
            //  in N).
            void moveOmpStep
            (
                std::vector<dosePathParticle*>& particles,
                std::vector<dosePathParticle::trackingData>& tds
            );


    public:


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dosePathCloud&) = delete;
};


} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
