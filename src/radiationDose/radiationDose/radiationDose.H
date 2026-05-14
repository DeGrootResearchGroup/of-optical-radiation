/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::functionObjects::radiationDose

Description
    Function object that integrates radiation dose along Lagrangian
    particle trajectories through a (typically frozen) flow field.

    Reads a velocity field U [m/s] and a fluence-rate field G [W/m^2],
    seeds particles through a configurable seeding model, and integrates
    their motion with optional stochastic turbulent dispersion through
    OpenFOAM's barycentric-tet particle tracker. Drift-free by
    construction: position is stored as barycentric coordinates within
    the current tet, never as Cartesian (x,y,z), so face crossings
    introduce no perpendicular floating-point error.

    Units convention:
    - Fluence rate input G is in W/m^2 (OpenFOAM SI standard).
    - Accumulated dose D is in mJ/cm^2 (UV reactor literature standard).
    - The unit conversion factor 0.1 (W/m^2 . s -> mJ/cm^2) is applied
      inside the integration kernel.
    - Inactivation rate kInact is in cm^2/mJ.

Usage
    \verbatim
    radiationDose
    {
        type            radiationDose;
        libs            ("libradiationDose.so");

        U               U;
        fluenceRate     G;
        seed            42;

        seeding
        {
            type        patchInjection;
            patches     (inlet);
            nParticles  10000;
        }

        dispersion
        {
            type        none;          // or discreteRandomWalk
        }

        termination
        {
            escapePatches   (outlet);
            maxTime         300;        // s; 0 disables
            maxDose         0;          // mJ/cm^2; 0 disables
            wallReflection  true;
        }

        integration
        {
            dtMax           0.01;       // s
            cflMax          0.5;
            maxOuterSteps   100000;     // safety cap on outer iterations
        }

        output
        {
            kInact          (0.1);      // cm^2/mJ for log-reduction reporting
            writeVtk        true;       // emit trajectories.vtk per execute()

            // Optional pre-write dose filter (inclusive bounds). A track
            // is written to trajectories.vtk only if its final dose D
            // satisfies vtkMinDose <= D <= vtkMaxDose. Either bound is
            // disabled by setting it to a negative value (the default).
            // CSV and summary are unaffected -- they always reflect the
            // full seeded population.
            vtkMinDose      -1;         // mJ/cm^2; <0 disables lower bound
            vtkMaxDose      -1;         // mJ/cm^2; <0 disables upper bound

            // Optional batching for the trajectory integration. Each
            // batch builds a fresh cloud of `batchSize` particles, runs
            // them to completion, writes one numbered VTK file
            // (trajectories_00001.vtk, ...), then drops the cloud
            // before the next batch starts. Peak memory scales as
            // O(batchSize) instead of O(nParticles). A
            // trajectories.pvd wrapper indexes all batch files so
            // ParaView can open the whole population as one dataset.
            // CSV and summary are aggregated across batches. Setting
            // batchSize to 0 (the default) or to >= nParticles disables
            // batching and writes a single unsuffixed trajectories.vtk
            // with no PVD wrapper (the v1 behaviour bit-for-bit).
            batchSize       0;          // particles per batch; 0 disables
        }
    }
    \endverbatim

SourceFiles
    radiationDose.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_radiationDose_H
#define functionObjects_radiationDose_H

#include "Switch.H"
#include "autoPtr.H"
#include "dispersionModel.H"
#include "dosePathCloud.H"
#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "motionModel.H"
#include "randomGenerator.H"
#include "seedingModel.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class radiationDose Declaration
\*---------------------------------------------------------------------------*/

class radiationDose
:
    public fvMeshFunctionObject
{
    // Forward declarations of writer-pipeline payloads. Defined fully
    // in the .C since both are implementation details:
    //   GatheredTracks   -- per-batch gather result (rank-major,
    //                       includes per-vertex data; lives only for
    //                       the duration of writing one batch's VTK).
    //   AggregatedTracks -- cross-batch accumulator (flat 1-D, summary
    //                       fields only, plus the per-batch VTK file
    //                       metadata needed by the PVD wrapper).
    struct GatheredTracks;
    struct AggregatedTracks;


    // Private Data — configuration

        word UName_;
        word GName_;
        label randomSeed_;

        autoPtr<dose::seedingModel> seeding_;

        //- Dispersion model dictionary (kept across read() so we can
        //  rebuild the model fresh per execute() — its internal
        //  state needs to start clean every run, but the
        //  configuration is read once).
        dictionary dispersionDict_;

        //- Motion model dictionary (same lifecycle as dispersionDict_).
        //  Absent in the user dict ⇒ default to `{ type tracer; }` so
        //  cases that don't care about particle inertia stay minimal.
        dictionary motionDict_;

        labelHashSet escapePatchIDs_;
        scalar maxTime_;
        scalar maxDose_;
        Switch wallReflection_;

        scalar dtMax_;
        scalar cflMax_;
        label maxOuterSteps_;

        scalarList kInact_;
        Switch writeVtk_;

        //- Optional inclusive bounds on per-track final dose D for
        //  VTK output. A negative value disables the corresponding
        //  bound (the default). Filter is applied at gather time --
        //  trajectories are still accumulated during integration
        //  because the final dose is only known after the run.
        scalar vtkMinDose_;
        scalar vtkMaxDose_;

        //- Particles per batch (0 disables batching). When enabled,
        //  the cloud is rebuilt and run once per ceil(nParticles /
        //  batchSize) batch, with one VTK file emitted per batch.
        label batchSize_;


    // Private Data — built per execute()

        //- Cloud of dose-tracking particles. Lives across execute()
        //  and write() so that the writer can iterate the cloud's
        //  particles for the per-track CSV.
        autoPtr<dose::dosePathCloud> cloud_;

        //- Cross-batch accumulator. Allocated on rank 0 only;
        //  populated incrementally inside execute() as each batch
        //  completes; consumed by write() to emit aggregated
        //  CSV/summary and the PVD wrapper.
        autoPtr<AggregatedTracks> aggregated_;


    // Private Member Functions

        //- Gather the local cloud's per-track and per-vertex data
        //  onto rank 0. In serial this is effectively a fill-and-
        //  return; in parallel, only the master rank's struct holds
        //  all ranks' data after the call. Vertex arrays are
        //  flattened across tracks (per rank) and the per-track
        //  nVerts list lets the consumer slice them.
        GatheredTracks gatherTracks() const;

        //- Append a freshly-gathered batch onto aggregated_. No-op
        //  on non-master ranks (aggregated_ is master-only).
        void appendBatch
        (
            const GatheredTracks& g,
            const word& vtkFileName
        );

        //- Run one batch: build cloud, seed with `seeds`, integrate,
        //  gather, write the per-batch VTK file, append to the
        //  cross-batch aggregator. Caller owns batch indexing and
        //  filename construction.
        void runBatch
        (
            const List<dose::trackPoint>& seeds,
            const word& vtkFileName,
            const fileName& outDir,
            const volVectorField& U,
            const volScalarField& G,
            randomGenerator& rng
        );

        void writeDoseCsv
        (
            const AggregatedTracks& a,
            const fileName& outDir
        ) const;
        void writeSummary
        (
            const AggregatedTracks& a,
            const fileName& outDir
        ) const;
        void writeVtkTrajectories
        (
            const GatheredTracks& g,
            const fileName& outFile
        ) const;
        void writePvdWrapper
        (
            const AggregatedTracks& a,
            const fileName& outDir
        ) const;


public:

    //- Runtime type information
    TypeName("radiationDose");


    // Constructors

        radiationDose
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        radiationDose(const radiationDose&) = delete;


    //- Destructor
    virtual ~radiationDose();


    // Member Functions

        virtual bool read(const dictionary&);
        virtual wordList fields() const;
        virtual bool execute();
        virtual bool write();


    // Member Operators

        void operator=(const radiationDose&) = delete;
};


} // End namespace functionObjects
} // End namespace Foam

#endif

// ************************************************************************* //
