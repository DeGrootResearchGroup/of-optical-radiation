/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dosePathCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(dosePathCloud, 0);
}
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::dose::dosePathCloud::dosePathCloud
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
)
:
    lagrangian::Cloud<dosePathParticle>(mesh, cloudName, false),
    mesh_(mesh),
    dtMax_(dtMax),
    cflMax_(cflMax),
    escapePatchIDs_(escapePatchIDs),
    maxTime_(maxTime),
    maxDose_(maxDose),
    wallReflection_(wallReflection),
    storeTrack_(storeTrack),
    dispersion_(std::move(dispersion))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dose::dosePathCloud::nActive() const
{
    label n = 0;
    forAllConstIter
    (
        typename Foam::lagrangian::Cloud<dosePathParticle>,
        *this,
        iter
    )
    {
        if (iter().active())
        {
            ++n;
        }
    }
    return n;
}


Foam::label Foam::dose::dosePathCloud::runToCompletion
(
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& GInterp,
    randomGenerator& rng,
    const label maxOuterSteps
)
{
    dosePathParticle::trackingData td(*this, UInterp, GInterp, rng);

    label step = 0;
    while (step < maxOuterSteps)
    {
        const label active = nActive();
        if (active == 0)
        {
            break;
        }
        if (debug && (step % 100 == 0))
        {
            Info<< "  outer step " << step
                << ": " << active << " active" << endl;
        }
        Cloud<dosePathParticle>::move(*this, td);
        ++step;
    }
    return step;
}


// ************************************************************************* //
