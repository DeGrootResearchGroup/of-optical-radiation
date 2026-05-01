/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::functionObjects::radiationDose

Description
    Function object that integrates radiation dose along Lagrangian
    particle trajectories through a (typically frozen) flow field.

    Reads a velocity field U [m/s] and a fluence-rate field G [W/m^2],
    seeds particles through a configurable seeding model, integrates
    their motion (with optional stochastic turbulent dispersion), and
    accumulates dose D along each trajectory.

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
        }

        integration
        {
            dtMax           0.01;       // s
            cflMax          0.5;
            maxStepsPerTrack 100000;
        }

        output
        {
            kInact          (0.1);      // cm^2/mJ for log-reduction reporting
        }
    }
    \endverbatim

SourceFiles
    radiationDose.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_radiationDose_H
#define functionObjects_radiationDose_H

#include "PtrList.H"
#include "Switch.H"
#include "autoPtr.H"
#include "dispersionModel.H"
#include "fvMeshFunctionObject.H"
#include "interpolation.H"
#include "HashSet.H"
#include "meshSearch.H"
#include "randomGenerator.H"
#include "seedingModel.H"
#include "track.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
    // Private Data — configuration

        //- Name of the velocity field [m/s]
        word UName_;

        //- Name of the fluence-rate field [W/m^2]
        word GName_;

        //- RNG seed for reproducible runs
        label randomSeed_;

        //- Seeding model (RTS-selected)
        autoPtr<dose::seedingModel> seeding_;

        //- Dispersion model (RTS-selected)
        autoPtr<dose::dispersionModel> dispersion_;

        //- Patches whose first hit terminates a track as "escaped"
        labelHashSet escapePatchIDs_;

        //- Soft termination: maximum integration time [s]; 0 disables
        scalar maxTime_;

        //- Soft termination: maximum accumulated dose [mJ/cm^2]; 0 disables
        scalar maxDose_;

        //- Specularly reflect off non-escape patches (default true).
        //  When false, a non-escape patch hit terminates the track as stuck.
        Switch wallReflection_;

        //- Maximum integration step [s]
        scalar dtMax_;

        //- CFL bound on per-step displacement
        scalar cflMax_;

        //- Hard cap on substeps per track (safety net)
        label maxStepsPerTrack_;

        //- Inactivation rate constants [cm^2/mJ] for log-reduction reporting
        scalarList kInact_;


    // Private Data — result of execute()

        //- Trajectories (one per seeded particle)
        PtrList<dose::track> tracks_;


    // Private Member Functions

        //- Integrate a single track to termination
        void integrateTrack
        (
            dose::track& tr,
            const interpolation<vector>& Uint,
            const interpolation<scalar>& Gint,
            const meshSearch& search,
            randomGenerator& rng
        ) const;

        //- Handle a step that crossed the domain boundary: locate the hit
        //  face, snap to it, classify endReason. When wallReflection_ is on
        //  and the hit patch is not in escapePatches, mirror the residual
        //  displacement about the face normal and return true so the outer
        //  loop continues. Returns false when the track has been terminated.
        bool handleBoundaryHit
        (
            dose::track& tr,
            const point& x_prev,
            label celli_prev,
            const vector& V,
            scalar dt,
            const interpolation<scalar>& Gint
        ) const;

        //- Write the per-track dose CSV
        void writeDoseCsv() const;

        //- Write the summary statistics file
        void writeSummary() const;


public:

    //- Conversion factor: G [W/m^2] * dt [s] -> D [mJ/cm^2]
    //  1 J/m^2 = 0.1 mJ/cm^2 (since 1 m^2 = 1e4 cm^2 and 1 J = 1e3 mJ)
    static constexpr scalar Wm2_s_to_mJcm2 = 0.1;


    //- Runtime type information
    TypeName("radiationDose");


    // Constructors

        //- Construct from Time and dictionary
        radiationDose
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        radiationDose(const radiationDose&) = delete;


    //- Destructor
    virtual ~radiationDose();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute (integrate trajectories and accumulate dose)
        virtual bool execute();

        //- Write outputs (dose CSV, summary)
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radiationDose&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
