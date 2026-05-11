/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Namespace
    Foam::dose

Description
    Sample a triple of independent N(0, 1) variates via the Box-Muller
    transform. Used by both DRW (resampling the eddy fluctuation u')
    and the inertial motion model (Brownian velocity kick).

    Three samples from four uniforms: the fourth uniform is wasted
    (Marsaglia-polar would avoid that, but Box-Muller is simpler and
    the cost is one extra log/sin/cos call per outer step, negligible
    against the integration work). r1 / r3 clamped to `small` so the
    log argument can never underflow to -inf when scalar01() returns
    exactly zero.

\*---------------------------------------------------------------------------*/

#ifndef dose_gaussianSample_H
#define dose_gaussianSample_H

#include "vector.H"
#include "randomGenerator.H"
#include "mathematicalConstants.H"
#include "scalar.H"

namespace Foam
{
namespace dose
{

inline vector gaussianTriple(randomGenerator& rng)
{
    const scalar twoPi = constant::mathematical::twoPi;

    const scalar r1 = max(rng.scalar01(), small);
    const scalar r2 = rng.scalar01();
    const scalar r3 = max(rng.scalar01(), small);
    const scalar r4 = rng.scalar01();

    const scalar mag1 = sqrt(-2.0*log(r1));
    const scalar mag3 = sqrt(-2.0*log(r3));

    return vector
    (
        mag1*cos(twoPi*r2),
        mag1*sin(twoPi*r2),
        mag3*cos(twoPi*r4)
    );
}

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
