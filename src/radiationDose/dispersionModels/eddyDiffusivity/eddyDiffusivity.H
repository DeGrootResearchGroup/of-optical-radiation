/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::eddyDiffusivity

Description
    The turbulence a random walk samples: the fluctuation scale
    sigma = sqrt(2k/3) per component, the eddy lifetime tau_e and the
    diffusivity K = sigma^2 tau_e / 2 = k tau_e / 3 of a walk of eddies
    held for tau_e, with its gradient. Shared by the dispersion models
    that need them, so each reads the same fields the same way.

    The dissipation comes from the turbulence model's own field: the
    rate epsilon of a k-epsilon model, giving tau_e = Cl k / epsilon, or
    the specific rate omega of a k-omega model, with epsilon = Cmu k omega
    and so tau_e = Cl / (Cmu omega). Name one of the two; epsilon is the
    default. tau_e is capped at tauEMax.

    Two ways to sample:

    - at(): k and the dissipation interpolated cell-point (linear in the
      particle's tetrahedron), sigma, tau_e and K formed from them, and
      grad(K) as the exact gradient of that K. Needs correct() to have
      built the interpolators from the registry's current fields; only
      reads them afterwards, so it is safe from concurrent threads.
    - cell(): k and the dissipation of one cell, as the uncorrected
      discrete random walk reads them.

Usage
    Entries read from the dispersion model's dictionary:

        k           k;            // optional, default "k"
        epsilon     epsilon;      // optional, default "epsilon"
        omega       omega;        // instead of epsilon; not both
        Cmu         0.09;         // optional, default 0.09 (with omega)
        Cl          0.15;         // optional, default 0.15
        tauEMax     100;          // optional, default 100 s

\*---------------------------------------------------------------------------*/

#ifndef dose_eddyDiffusivity_H
#define dose_eddyDiffusivity_H

#include "autoPtr.H"
#include "barycentric.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "interpolationCellPoint.H"
#include "tetIndices.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class eddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class eddyDiffusivity
{
public:

    //- The turbulence at a point
    struct sample
    {
        //- Per-component fluctuation scale sqrt(2k/3) [m/s]
        scalar sigma;

        //- Eddy lifetime [s]
        scalar tau;

        //- Gradient of the eddy lifetime [s/m]
        vector gradTau;

        //- Diffusivity K = k tau_e / 3 [m^2/s]
        scalar K;

        //- Gradient of K [m/s]
        vector gradK;

        //- Gradient of sigma [1/s]: grad(k) / (3 sigma), zero where k is
        //  nil
        vector gradSigma;
    };


private:

    // Private Data

        const fvMesh& mesh_;

        //- Name of the turbulent kinetic energy field [m^2/s^2]
        word kName_;

        //- True when the dissipation field is omega rather than epsilon
        bool omega_;

        //- Name of the dissipation field: epsilon [m^2/s^3] or omega [1/s]
        word dissipationName_;

        //- epsilon = Cmu * k * omega (used only with omega)
        scalar Cmu_;

        //- Eddy lifetime constant: tau_e = Cl * k / epsilon
        scalar Cl_;

        //- Hard upper bound on the eddy lifetime [s]. In stagnant
        //  zones (eps -> 0 with finite k) the raw `Cl k/eps` is
        //  unbounded, locking the same fluctuation `u'` onto a
        //  particle for the entire trajectory; capping at `tauEMax`
        //  forces a resample on a sensible timescale. Default is a
        //  generous 100 s: large enough to leave normal flow eddies
        //  alone (a UV reactor's largest eddies are O(seconds)) but
        //  finite enough to bound the stagnant degenerate case.
        scalar tauEMax_;

        //- Linear-in-tetrahedron interpolators, built by correct()
        autoPtr<interpolationCellPoint<scalar>> kInterp_;
        autoPtr<interpolationCellPoint<scalar>> dissInterp_;


public:

    // Constructors

        eddyDiffusivity(const dictionary& dict, const fvMesh& mesh);


    // Member Functions

        //- The dissipation field the dict names: its "omega" entry if
        //  it has one, else its "epsilon" entry (default "epsilon").
        //  Naming both is an error.
        static word dissipationName(const dictionary& dict);

        //- The fields the dict's entries name: k and the dissipation
        static wordList requiredFields(const dictionary& dict);

        wordList requiredFields() const
        {
            return wordList({kName_, dissipationName_});
        }

        const word& dissipationName() const
        {
            return dissipationName_;
        }

        //- Build the interpolators from the registry's current fields
        void correct();

        //- Interpolated turbulence at the coordinates in the tetrahedron
        sample at
        (
            const barycentric& coordinates,
            const tetIndices& tetIs
        ) const;

        //- One cell's k [m^2/s^2], clipped at zero, and eddy lifetime [s]
        void cell(label celli, scalar& k, scalar& tau) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
