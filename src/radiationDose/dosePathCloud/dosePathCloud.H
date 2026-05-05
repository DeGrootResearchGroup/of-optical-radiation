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
            inline dispersionModel& dispersion() { return *dispersion_; }
            inline const dispersionModel& dispersion() const
            {
                return *dispersion_;
            }


        // Edit

            //- Drive Cloud::move() in a loop until every particle is
            //  no longer active (escaped / timedOut / stuck /
            //  terminated) or maxOuterSteps is reached as a safety
            //  net. Returns the number of outer steps actually run.
            label runToCompletion
            (
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& GInterp,
                randomGenerator& rng,
                const label maxOuterSteps
            );

            //- Number of particles still in the active state
            //  (endReason == active). Counted via a forAllConstIter.
            label nActive() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dosePathCloud&) = delete;
};


} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
