/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

// Inline constructor for the trackingData. Defined here rather than in
// the class body to break the include cycle with dosePathCloud.

inline Foam::dose::dosePathParticle::trackingData::trackingData
(
    const dosePathCloud& cloud,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& GInterp,
    randomGenerator& rng
)
:
    particle::trackingData(cloud),
    UInterp_(UInterp),
    GInterp_(GInterp),
    rng_(rng)
{}


// ************************************************************************* //
