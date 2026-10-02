/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
\*---------------------------------------------------------------------------*/

#include "eddyDiffusivity.H"
#include "triFace.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * * //

namespace
{

// Gradient of a cell-point interpolated field in one tetrahedron. The
// interpolation is linear there, with the cell value at the cell centre
// (vertex a) and point values at the three face points (b, c, d), so
// the gradient is constant: with e1 = b - a, e2 = c - a, e3 = d - a and
// the value differences dv_i,
//     grad = (dv_1 e2 x e3 + dv_2 e3 x e1 + dv_3 e1 x e2) / (e1 . e2 x e3).
// A degenerate tetrahedron has no gradient; zero is returned.
Foam::vector tetGradient
(
    const Foam::interpolationCellPoint<Foam::scalar>& field,
    const Foam::polyMesh& mesh,
    const Foam::tetIndices& tetIs
)
{
    using namespace Foam;

    const triFace tri = tetIs.faceTriIs(mesh);
    const tetPointRef t = tetIs.tet(mesh);

    const vector e1 = t.b() - t.a();
    const vector e2 = t.c() - t.a();
    const vector e3 = t.d() - t.a();
    const vector n1 = e2 ^ e3;
    const scalar det = e1 & n1;

    if (mag(det) <= vSmall)
    {
        return vector::zero;
    }

    const scalar va = field.psi()[tetIs.cell()];
    const scalarField& vp = field.psip();

    return
    (
        (vp[tri[0]] - va)*n1
      + (vp[tri[1]] - va)*(e3 ^ e1)
      + (vp[tri[2]] - va)*(e1 ^ e2)
    )/det;
}

} // End anonymous namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::dose::eddyDiffusivity::dissipationName
(
    const dictionary& dict
)
{
    if (dict.found("omega") && dict.found("epsilon"))
    {
        FatalIOErrorInFunction(dict)
            << "name either epsilon or omega, not both"
            << exit(FatalIOError);
    }
    return dict.found("omega")
        ? dict.lookup<word>("omega")
        : dict.lookupOrDefault<word>("epsilon", "epsilon");
}


Foam::wordList Foam::dose::eddyDiffusivity::requiredFields
(
    const dictionary& dict
)
{
    return wordList
    ({
        dict.lookupOrDefault<word>("k", "k"),
        dissipationName(dict)
    });
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dose::eddyDiffusivity::eddyDiffusivity
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    omega_(dict.found("omega")),
    dissipationName_(dissipationName(dict)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Cl_(dict.lookupOrDefault<scalar>("Cl", 0.15)),
    tauEMax_(dict.lookupOrDefault<scalar>("tauEMax", 100.0)),
    kInterp_(),
    dissInterp_()
{
    if (Cmu_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Cmu must be > 0, got " << Cmu_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dose::eddyDiffusivity::correct()
{
    kInterp_.reset
    (
        new interpolationCellPoint<scalar>
        (
            mesh_.lookupObject<volScalarField>(kName_)
        )
    );
    dissInterp_.reset
    (
        new interpolationCellPoint<scalar>
        (
            mesh_.lookupObject<volScalarField>(dissipationName_)
        )
    );
}


Foam::dose::eddyDiffusivity::sample Foam::dose::eddyDiffusivity::at
(
    const barycentric& coordinates,
    const tetIndices& tetIs
) const
{
    if (!kInterp_.valid() || !dissInterp_.valid())
    {
        FatalErrorInFunction
            << "correct() has not been called, so there are no "
            << kName_ << " / " << dissipationName_
            << " interpolators to sample" << exit(FatalError);
    }

    scalar k = kInterp_->interpolate(coordinates, tetIs);
    vector gradk = tetGradient(kInterp_(), mesh_, tetIs);
    if (k <= 0)
    {
        k = 0;
        gradk = vector::zero;
    }

    scalar diss = dissInterp_->interpolate(coordinates, tetIs);
    vector gradDiss = tetGradient(dissInterp_(), mesh_, tetIs);
    if (diss <= small)
    {
        diss = small;
        gradDiss = vector::zero;
    }

    // tau_e and its gradient, from Cl / (Cmu omega) or Cl k / epsilon
    scalar tau;
    vector gradTau;
    if (omega_)
    {
        tau = Cl_/(Cmu_*diss);
        gradTau = -(tau/diss)*gradDiss;
    }
    else
    {
        tau = Cl_*k/diss;
        gradTau = (Cl_/diss)*(gradk - (k/diss)*gradDiss);
    }
    if (tau >= tauEMax_)
    {
        tau = tauEMax_;
        gradTau = vector::zero;
    }

    // K = sigma^2 tau / 2 = k tau / 3
    sample s;
    s.sigma = sqrt(2.0/3.0*k);
    s.tau = tau;
    s.gradTau = gradTau;
    s.K = k*tau/3.0;
    s.gradK = (tau*gradk + k*gradTau)/3.0;
    s.gradSigma = s.sigma > 0 ? gradk/(3.0*s.sigma) : vector::zero;
    return s;
}


void Foam::dose::eddyDiffusivity::cell
(
    label celli,
    scalar& k,
    scalar& tau
) const
{
    const volScalarField& kField =
        mesh_.lookupObject<volScalarField>(kName_);
    const volScalarField& dissField =
        mesh_.lookupObject<volScalarField>(dissipationName_);

    k = max(kField[celli], scalar(0));
    const scalar diss = max(dissField[celli], small);

    // Cl k / epsilon, which with epsilon = Cmu k omega is Cl / (Cmu omega)
    tau = min(omega_ ? Cl_/(Cmu_*diss) : Cl_*k/diss, tauEMax_);
}


// ************************************************************************* //
