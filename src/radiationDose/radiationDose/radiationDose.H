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
    // Forward declaration of the per-rank gathered track payload.
    // Defined fully in the .C since it's an implementation detail of
    // the writer pipeline.
    struct GatheredTracks;


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


    // Private Data — built per execute()

        //- Cloud of dose-tracking particles. Lives across execute()
        //  and write() so that the writer can iterate the cloud's
        //  particles for the per-track CSV.
        autoPtr<dose::dosePathCloud> cloud_;


    // Private Member Functions

        //- Gather the local cloud's per-track and per-vertex data
        //  onto rank 0. In serial this is effectively a fill-and-
        //  return; in parallel, only the master rank's struct holds
        //  all ranks' data after the call. Vertex arrays are
        //  flattened across tracks (per rank) and the per-track
        //  nVerts list lets the consumer slice them.
        GatheredTracks gatherTracks() const;

        void writeDoseCsv
        (
            const GatheredTracks& g,
            const fileName& outDir
        ) const;
        void writeSummary
        (
            const GatheredTracks& g,
            const fileName& outDir
        ) const;
        void writeVtkTrajectories
        (
            const GatheredTracks& g,
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
