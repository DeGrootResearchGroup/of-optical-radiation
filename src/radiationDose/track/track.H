/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::track

Description
    A single particle trajectory: a sequence of trackPoint vertices plus
    metadata (track ID, termination reason).

\*---------------------------------------------------------------------------*/

#ifndef dose_track_H
#define dose_track_H

#include "DynamicList.H"
#include "NamedEnum.H"
#include "autoPtr.H"
#include "dispersionModel.H"
#include "trackPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                            Class track Declaration
\*---------------------------------------------------------------------------*/

class track
{
public:

    //- Reason a track stopped accumulating dose
    enum class endReason
    {
        active,        // still tracking (initial / mid-run)
        escaped,       // hit a configured escape patch
        timedOut,      // exceeded max integration time
        stuck,         // step size collapsed (degenerate cell or zero U)
        leftDomain,    // walked off mesh without a patch hit (parallel/error)
        terminated     // generic termination
    };

    //- Names for endReason values (for diagnostics / output)
    static const NamedEnum<endReason, 6> endReasonNames;


private:

    // Private Data

        //- Track ID (unique within a single execute() invocation)
        label id_;

        //- Vertices along the trajectory
        DynamicList<trackPoint> points_;

        //- Why the track stopped accumulating
        endReason end_;

        //- Per-track dispersion-model state. Owned by the track so each
        //  thread can mutate its own particles' state without locking.
        //  Set by the function object at seed time via the dispersion
        //  model's newState() factory.
        autoPtr<dispersionModel::State> dispState_;


public:

    // Constructors

        //- Construct empty with the given ID
        explicit track(label id, label initialCapacity = 64)
        :
            id_(id),
            points_(initialCapacity),
            end_(endReason::active),
            dispState_()
        {}


    // Member Functions

        // Access

            inline label id() const { return id_; }
            inline const DynamicList<trackPoint>& points() const { return points_; }
            inline label size() const { return points_.size(); }
            inline endReason end() const { return end_; }

            //- The most recently appended point (UB if empty)
            inline const trackPoint& back() const { return points_.last(); }

            //- Final accumulated dose [mJ/cm^2] (0 if empty)
            inline scalar finalDose() const
            {
                return points_.empty() ? 0 : points_.last().D;
            }

            //- Total residence time [s] (0 if empty)
            inline scalar finalTime() const
            {
                return points_.empty() ? 0 : points_.last().t;
            }

            //- Mutable per-track dispersion state. Caller owns thread
            //  exclusivity (one track per thread at a time).
            inline dispersionModel::State& dispState()
            {
                return *dispState_;
            }

        // Edit

            //- Append a vertex
            inline void append(const trackPoint& p) { points_.append(p); }

            //- Mark the termination reason
            inline void setEnd(endReason r) { end_ = r; }

            //- Install a per-track dispersion state (called once at
            //  seed time, before the integration loop).
            inline void setDispState(autoPtr<dispersionModel::State> s)
            {
                dispState_ = std::move(s);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
