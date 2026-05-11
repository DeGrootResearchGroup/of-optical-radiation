/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::perTrackState

Description
    Empty polymorphic base for "per-track state" slots that
    dosePathParticle carries on behalf of its installed motion and
    dispersion models. A concrete model defines its own derived state
    (e.g. discreteRandomWalk::DRWState carries the current eddy
    fluctuation and remaining lifetime); the particle holds one
    autoPtr<perTrackState> per slot, and the model downcasts on
    access. Models that need no per-track memory use this base
    directly.

    motionModel::State and dispersionModel::State are aliases of
    this type. They used to be near-identical empty classes nested
    inside their respective bases; sharing one base eliminates the
    duplication without changing any user-facing API.

\*---------------------------------------------------------------------------*/

#ifndef dose_perTrackState_H
#define dose_perTrackState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

class perTrackState
{
public:
    virtual ~perTrackState() = default;
};

} // End namespace dose
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
