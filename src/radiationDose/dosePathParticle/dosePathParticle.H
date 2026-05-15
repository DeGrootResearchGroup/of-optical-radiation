/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::dosePathParticle

Description
    A particle that walks through a frozen flow field and accumulates
    radiation dose D = integral(G dt) along its trajectory.

    Derives from Foam::particle so that position is stored in
    barycentric tet coordinates rather than Cartesian (x,y,z). This
    eliminates floating-point drift during long trajectories: after a
    tet-face crossing, one barycentric coordinate is exactly zero, so
    there is no perpendicular numerical error to compound. Cartesian
    position is reconstructed on demand from the barycentric coordinates
    and the current tet's vertices.

    Each call to move() advances the particle by one outer time step of
    duration dtMax (CFL-bounded), sampling the flow velocity once at
    the start of the step and adding a turbulent fluctuation drawn from
    the cloud's dispersion model. Within the outer step the particle
    traverses an arbitrary number of cell faces via OpenFOAM's
    trackToFace machinery, accruing dose at each segment.

    Wall reflection: hitWallPatch() specularly reflects V_ about the
    wall normal. The particle stays in the same cell after reflection
    and continues the outer step with the new velocity.

    Escape: a hit on a patch listed in the cloud's escapePatchIDs marks
    endReason_ = escaped and stops further tracking.

\*---------------------------------------------------------------------------*/

#ifndef dose_dosePathParticle_H
#define dose_dosePathParticle_H

#include "Cloud.H"
#include "DynamicList.H"
#include "NamedEnum.H"
#include "autoPtr.H"
#include "contiguous.H"
#include "dispersionModel.H"
#include "interpolationCellPoint.H"
#include "motionModel.H"
#include "particle.H"
#include "trackPoint.H"

namespace Foam
{

namespace dose
{

class dosePathCloud;

class dosePathParticle;

Ostream& operator<<(Ostream&, const dosePathParticle&);


/*---------------------------------------------------------------------------*\
                       Class dosePathParticle Declaration
\*---------------------------------------------------------------------------*/

class dosePathParticle
:
    public particle
{
public:

    //- Why a track stopped accumulating dose. Mirrors the v0.2 enum
    //  but without the speculative "stalled" / "drifted" classes —
    //  with barycentric tracking, drift is impossible by construction
    //  and stall is a downstream concern (a recirculation cell is
    //  caught by maxTime, which is the physically meaningful bound).
    enum class endReason
    {
        active,        // still tracking
        escaped,       // hit a configured escape patch
        timedOut,      // exceeded cloud.maxTime()
        stuck,         // step size collapsed (zero V), or non-escape
                       // boundary hit when wallReflection is disabled
        terminated     // exceeded cloud.maxDose()
    };

    static const NamedEnum<endReason, 5> endReasonNames;


private:

    // Private Data — bytewise-trivial members must come first for
    // contiguous-IO support (see contiguous<dosePathParticle> below).

        //- Particle velocity [m/s] at the END of the most recent outer
        //  step (the OU-update result for inertial; equal to U + u'
        //  for tracer). Carried across outer steps as V_old for the
        //  next motion-model advance() call.
        vector V_;

        //- Mean velocity [m/s] over the current outer step (Δx/Δt).
        //  This is the vector used by the inner loop for
        //  trackToAndHitFace, and the one reflected by hitWallPatch
        //  during inner-step wall interactions. For tracer V_disp_ ==
        //  V_; for inertial they differ by the OU phi factor.
        vector V_disp_;

        //- Accumulated dose [mJ/cm^2]
        scalar D_;

        //- Integration time elapsed [s]
        scalar t_;

        //- Termination reason (initial: active)
        endReason endReason_;


    // Non-trivial members — not contiguous, not part of the on-disk
    // payload (the dose tracker is purely in-memory).

        //- Per-particle dispersion state. Owned by the particle so each
        //  thread / processor can mutate its own particles' state
        //  without coordination.
        autoPtr<dispersionModel::State> dispState_;

        //- Per-particle motion-model state (mirrors dispState_).
        //  Empty for tracer / drag-only inertial; reserved for
        //  motion models that need lagged per-particle quantities.
        autoPtr<motionModel::State> motionState_;

        //- Trajectory: every position the particle has visited at the
        //  end of an outer step. Controlled by the cloud's storeTrack
        //  flag. Transmitted across processor patches via the particle
        //  (de)serialisation operators so the VTK writer sees the
        //  full pre-handoff trajectory after parallel transfer.
        DynamicList<trackPoint> points_;


public:

    friend class Foam::lagrangian::Cloud<dosePathParticle>;


    //- Class used to pass tracking data to trackToAndHitFace.
    //  Carries references to the field interpolators, the dispersion
    //  model state container, and the shared RNG.
    class trackingData
    :
        public particle::trackingData
    {
        const interpolationCellPoint<vector>& UInterp_;
        const interpolationCellPoint<scalar>& GInterp_;
        randomGenerator& rng_;

    public:

        trackingData
        (
            const dosePathCloud& cloud,
            const interpolationCellPoint<vector>& UInterp,
            const interpolationCellPoint<scalar>& GInterp,
            randomGenerator& rng
        );

        const interpolationCellPoint<vector>& UInterp() const
        {
            return UInterp_;
        }
        const interpolationCellPoint<scalar>& GInterp() const
        {
            return GInterp_;
        }
        randomGenerator& rng() const { return rng_; }
    };


    //- Conversion factor: G [W/m^2] * dt [s] -> D [mJ/cm^2]
    //  1 J/m^2 = 0.1 mJ/cm^2 (since 1 m^2 = 1e4 cm^2 and 1 J = 1e3 mJ)
    static constexpr scalar Wm2_s_to_mJcm2 = 0.1;


    //- Runtime type information
    TypeName("dosePathParticle");


    // Constructors

        //- Construct from a position and a cell, locating the
        //  surrounding tet via the meshSearch engine.
        dosePathParticle
        (
            const meshSearch& searchEngine,
            const vector& position,
            const label celli,
            label& nLocateBoundaryHits
        );

        //- Construct from Istream. The dose tracker does not persist
        //  particles to disk in v0.3, but Cloud<...> requires this
        //  constructor to instantiate; we read just the base particle
        //  state and zero our own fields.
        dosePathParticle(Istream& is, bool readFields = true);

        //- Copy constructor
        dosePathParticle(const dosePathParticle& p);

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new dosePathParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<dosePathParticle> New(Istream& is)
        {
            return autoPtr<dosePathParticle>(new dosePathParticle(is));
        }


    // Member Functions

        // Access

            inline const vector& V() const { return V_; }
            inline scalar D() const { return D_; }
            inline scalar t() const { return t_; }
            inline endReason end() const { return endReason_; }
            inline bool active() const { return endReason_ == endReason::active; }
            inline const DynamicList<trackPoint>& points() const
            {
                return points_;
            }


        // Edit

            //- Install a per-particle dispersion state (called once at
            //  seed time, before the first move() invocation).
            inline void setDispState(autoPtr<dispersionModel::State> s)
            {
                dispState_ = std::move(s);
            }

            //- Install a per-particle motion-model state (called once
            //  at seed time).
            inline void setMotionState(autoPtr<motionModel::State> s)
            {
                motionState_ = std::move(s);
            }

            //- Append a track point. Called by the cloud at seed time
            //  to record the initial position, and (optionally) by
            //  move() at the end of each outer step.
            inline void appendPoint(const trackPoint& p)
            {
                points_.append(p);
            }


        // Tracking

            //- Advance the particle by one outer step (dtMax of physical
            //  time, CFL-bounded). Returns true if the particle should
            //  be kept in the cloud, false if deleted.
            bool move(dosePathCloud& cloud, trackingData& td);


        // Patch interactions

            //- Specular reflection of V_ about the wall normal.
            //  When the cloud's wallReflection flag is false, this
            //  marks the particle as stuck instead.
            void hitWallPatch(dosePathCloud& cloud, trackingData& td);

            //- Default boundary-patch interaction. Used for any patch
            //  that isn't specifically a wall, processor, cyclic,
            //  symmetry, etc. — including the user's escape patches.
            //  Sets endReason_ to escaped (if listed) or stuck
            //  (otherwise).
            void hitBasicPatch(dosePathCloud& cloud, trackingData& td);


    // I-O

        //- Read fields from disk into a cloud. The dose tracker is
        //  purely in-memory; this is a no-op stub required by Cloud<>.
        static void readFields
        (
            Foam::lagrangian::Cloud<dosePathParticle>& c
        );

        //- Write fields. Same: stub. Per-track output goes through
        //  the function object's CSV writer, not the Cloud I/O path.
        static void writeFields
        (
            const Foam::lagrangian::Cloud<dosePathParticle>& c
        );


    // Friend Operators

        friend Ostream& operator<<(Ostream&, const dosePathParticle&);
};


} // End namespace dose


// The contiguous trait must be in the Foam namespace (where the
// primary template lives). The dose tracker has non-trivial members
// (DynamicList<trackPoint>, autoPtr<dispersionModel::State>), so
// conservative answer is false; the slow generic Pstream path
// handles serialisation.
template<>
inline bool contiguous<dose::dosePathParticle>()
{
    return false;
}


} // End namespace Foam


#include "dosePathParticleI.H"


#endif

// ************************************************************************* //
