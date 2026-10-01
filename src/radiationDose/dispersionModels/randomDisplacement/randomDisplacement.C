/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "randomDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "gaussianSample.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{
    defineTypeNameAndDebug(randomDisplacement, 0);
    addToRunTimeSelectionTable(dispersionModel, randomDisplacement, dictionary);

    namespace
    {
        const dispersionModel::addRequiredFields _rdReqFields
        (
            "randomDisplacement",
            &randomDisplacement::requiredFields
        );
    }
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::randomDisplacement::randomDisplacement
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dispersionModel(dict, mesh),
    turbulence_(dict, mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dose::randomDisplacement::correct()
{
    turbulence_.correct();
}


void Foam::dose::randomDisplacement::beginStep
(
    dispersionModel::State& state,
    const barycentric& coordinates,
    const tetIndices& tetIs,
    scalar dtMax,
    randomGenerator& rng
) const
{
    RDState& s = dynamic_cast<RDState&>(state);

    const eddyDiffusivity::sample t = turbulence_.at(coordinates, tetIs);
    s.K = t.K;
    s.drift = t.gradK;
    s.xi = gaussianTriple(rng);
}


Foam::vector Foam::dose::randomDisplacement::fluctuation
(
    const dispersionModel::State& state,
    scalar dt
) const
{
    const RDState& s = dynamic_cast<const RDState&>(state);
    return s.drift + sqrt(2*s.K/dt)*s.xi;
}


Foam::vector Foam::dose::randomDisplacement::cflFluctuation
(
    const dispersionModel::State& state,
    scalar dt
) const
{
    return dynamic_cast<const RDState&>(state).drift;
}


Foam::scalar Foam::dose::randomDisplacement::maxStep
(
    const dispersionModel::State& state,
    scalar maxDisplacement
) const
{
    const RDState& s = dynamic_cast<const RDState&>(state);
    return s.K > vSmall ? sqr(maxDisplacement)/(2*s.K) : great;
}


// ************************************************************************* //
