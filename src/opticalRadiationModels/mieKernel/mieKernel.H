/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::optical::mieKernel

Description
    Mie scattering kernel for an isotropic homogeneous sphere of arbitrary
    size parameter and complex refractive index. Implements the
    Bohren-Huffman BHMIE recurrence: downward recurrence on the
    logarithmic derivative D_n(mx) = psi_n'(mx)/psi_n(mx), upward
    recurrence on the Riccati-Bessel functions psi_n(x), chi_n(x), and
    closed-form evaluation of the Mie coefficients a_n, b_n. The series
    is truncated at

        N_max = ceil(x + 4 x^(1/3) + 2)

    (Wiscombe 1980), with the downward recurrence started at
    N_start = max(N_max, |mx|) + 15 to ensure convergence to the correct
    physical D_n.

    Conventions:
      - x  = 2 pi r n_medium / lambda_vacuum         (size parameter)
      - m  = m_particle / n_medium                    (relative index, complex)
        Imaginary part of m_particle is non-negative for absorbing
        material in the exp(-i omega t) time convention used here.
      - Phase function returns |S_1(mu)|^2 + |S_2(mu)|^2 (proportional
        to dC_sca/dOmega for unpolarised incident light); absolute
        scaling is the caller's responsibility. The discrete-ordinates
        in-scatter pipeline (HenyeyGreensteinModel, rayleighModel)
        renormalises so the per-row sum is unity, so only the *shape*
        of the phase function matters in the table-construction loop.

References
    [1] Bohren & Huffman, "Absorption and Scattering of Light by Small
        Particles", Wiley, 1983, Appendix A (BHMIE).
    [2] Wiscombe, W. J., "Improved Mie scattering algorithms",
        Appl. Opt. 19, 1505-1509 (1980).

SourceFiles
    mieKernel.C

\*---------------------------------------------------------------------------*/

#ifndef optical_mieKernel_H
#define optical_mieKernel_H

#include "scalar.H"
#include "label.H"
#include "List.H"

#include <complex>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                          Class mieKernel Declaration
\*---------------------------------------------------------------------------*/

class mieKernel
{
public:

    typedef std::complex<scalar> complex;

private:

    // Private data

        //- Size parameter
        scalar x_;

        //- Relative complex refractive index
        complex m_;

        //- Number of terms retained in the Mie series
        label N_;

        //- Mie coefficients a_n, b_n  (size N_, indexed n = 1..N_)
        List<complex> an_;
        List<complex> bn_;

        //- Extinction efficiency (cross-section / pi r^2)
        scalar Qext_;

        //- Scattering efficiency
        scalar Qsca_;

        //- Asymmetry parameter <cos theta>
        scalar g_;

public:

    // Constructors

        //- Construct from size parameter and relative complex index. The
        //  Mie series is evaluated immediately; access methods are O(1)
        //  except phaseIntensity() which is O(N).
        mieKernel(const scalar x, const complex& m);


    // Member Functions

        //- Number of Mie terms retained (Wiscombe truncation)
        inline label N() const { return N_; }

        //- Extinction efficiency
        inline scalar Qext() const { return Qext_; }

        //- Scattering efficiency
        inline scalar Qsca() const { return Qsca_; }

        //- Absorption efficiency
        inline scalar Qabs() const { return Qext_ - Qsca_; }

        //- Asymmetry parameter g = <cos theta>
        inline scalar g() const { return g_; }

        //- Returns |S_1(mu)|^2 + |S_2(mu)|^2 for cosine of scattering
        //  angle mu = cos(theta). Proportional to dC_sca/dOmega for
        //  unpolarised incident light. Absolute normalisation is the
        //  caller's responsibility.
        scalar phaseIntensity(const scalar mu) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
