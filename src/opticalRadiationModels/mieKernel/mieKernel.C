/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "mieKernel.H"
#include "error.H"

#include <algorithm>
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::mieKernel::mieKernel
(
    const scalar x,
    const complex& m
)
:
    x_(x),
    m_(m),
    N_(0),
    an_(),
    bn_(),
    Qext_(0),
    Qsca_(0),
    g_(0)
{
    if (x_ <= 0)
    {
        FatalErrorIn("mieKernel::mieKernel(scalar, complex)")
            << "size parameter x must be positive (got " << x_ << ")"
            << exit(FatalError);
    }
    if (std::abs(m_) <= 0)
    {
        FatalErrorIn("mieKernel::mieKernel(scalar, complex)")
            << "relative refractive index m must be non-zero"
            << exit(FatalError);
    }

    // Wiscombe (1980) truncation: N_max = ceil(x + 4 x^(1/3) + 2). At
    // small x this is dominated by the +2 term; for x = 0.1 it gives
    // N_max = 5, more than enough for the dipole limit.
    N_ = label(std::ceil(x_ + 4.0*std::cbrt(x_) + 2.0));

    const complex y = m_*x_;

    // Downward-recurrence start point. Bohren-Huffman BHMIE uses
    //     N_start = max(N_max, |y|) + 15
    // which gives ~15 buffer iterations for the recurrence to converge
    // to the correct physical D_n from the seed D_{N_start} = 0.
    const label Nstart =
        label(std::max(scalar(N_), scalar(std::abs(y)))) + 15;

    // Downward recurrence for the logarithmic derivative
    //   D_{n-1}(y) = n/y - 1/(D_n(y) + n/y)
    // Indexed 0..Nstart; only entries 1..N_ are used afterwards.
    List<complex> D(Nstart + 1, complex(0, 0));
    for (label n = Nstart; n >= 1; --n)
    {
        const complex rn = complex(scalar(n), 0)/y;
        D[n - 1] = rn - complex(1, 0)/(D[n] + rn);
    }

    // Sized N_+1 so we can index 1..N_ directly; index 0 unused.
    an_.setSize(N_ + 1);
    bn_.setSize(N_ + 1);
    an_[0] = bn_[0] = complex(0, 0);

    // Upward recurrence for the real Riccati-Bessel functions
    //   psi_n(x) and chi_n(x).
    // Initial values:
    //   psi_{-1}(x) = cos(x), psi_0(x) = sin(x)
    //   chi_{-1}(x) = -sin(x), chi_0(x) =  cos(x)
    // Recurrence (same for both):
    //   f_n = (2n - 1)/x * f_{n-1} - f_{n-2}
    // The Riccati-Hankel function is xi_n = psi_n - i chi_n.
    scalar psi_nm2 = std::cos(x_);
    scalar psi_nm1 = std::sin(x_);
    scalar chi_nm2 = -std::sin(x_);
    scalar chi_nm1 = std::cos(x_);

    for (label n = 1; n <= N_; ++n)
    {
        const scalar fn = scalar(n);

        const scalar psi_n = (2.0*fn - 1.0)/x_ * psi_nm1 - psi_nm2;
        const scalar chi_n = (2.0*fn - 1.0)/x_ * chi_nm1 - chi_nm2;

        const complex xi_n  (psi_n,   -chi_n);
        const complex xi_nm1(psi_nm1, -chi_nm1);

        const complex Dn = D[n];

        const complex alpha = Dn/m_ + complex(fn/x_, 0);
        const complex beta  = m_*Dn + complex(fn/x_, 0);

        an_[n] =
            (alpha*psi_n - psi_nm1)
           /(alpha*xi_n  - xi_nm1);
        bn_[n] =
            (beta*psi_n - psi_nm1)
           /(beta*xi_n  - xi_nm1);

        const scalar twoNplus1 = 2.0*fn + 1.0;
        Qext_ += twoNplus1*(an_[n].real() + bn_[n].real());
        Qsca_ += twoNplus1*(std::norm(an_[n]) + std::norm(bn_[n]));

        psi_nm2 = psi_nm1;
        psi_nm1 = psi_n;
        chi_nm2 = chi_nm1;
        chi_nm1 = chi_n;
    }

    Qext_ *= 2.0/(x_*x_);
    Qsca_ *= 2.0/(x_*x_);

    // Asymmetry parameter g = <cos theta>:
    //   x^2 Q_sca g = 4 sum_{n=1}^{N-1} n(n+2)/(n+1) Re(a_n a*_{n+1} + b_n b*_{n+1})
    //               + 4 sum_{n=1}^{N}   (2n+1)/(n(n+1)) Re(a_n b*_n)
    scalar gNum = 0;
    for (label n = 1; n < N_; ++n)
    {
        const scalar fn = scalar(n);
        const scalar w = fn*(fn + 2.0)/(fn + 1.0);
        gNum += w*
        (
            (an_[n]*std::conj(an_[n + 1])).real()
          + (bn_[n]*std::conj(bn_[n + 1])).real()
        );
    }
    for (label n = 1; n <= N_; ++n)
    {
        const scalar fn = scalar(n);
        const scalar w = (2.0*fn + 1.0)/(fn*(fn + 1.0));
        gNum += w*(an_[n]*std::conj(bn_[n])).real();
    }
    if (Qsca_ > 0)
    {
        g_ = 4.0/(x_*x_*Qsca_)*gNum;
    }
    else
    {
        g_ = 0;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::optical::mieKernel::phaseIntensity(const scalar mu) const
{
    // Angular functions pi_n(mu), tau_n(mu) by upward recurrence.
    //   pi_1 = 1
    //   pi_n = (2n-1)/(n-1) mu pi_{n-1} - n/(n-1) pi_{n-2}    (n >= 2)
    //   tau_n = n mu pi_n - (n+1) pi_{n-1}
    //
    // S_1(mu) = sum_{n=1}^{N} (2n+1)/(n(n+1)) (a_n pi_n + b_n tau_n)
    // S_2(mu) = sum_{n=1}^{N} (2n+1)/(n(n+1)) (a_n tau_n + b_n pi_n)
    //
    // Returns |S_1|^2 + |S_2|^2 (proportional to dC_sca/dOmega for
    // unpolarised incident light).
    scalar pi_nm1 = 0;   // pi_0
    scalar pi_n   = 1;   // pi_1
    complex S1(0, 0);
    complex S2(0, 0);

    for (label n = 1; n <= N_; ++n)
    {
        const scalar fn = scalar(n);
        const scalar tau_n = fn*mu*pi_n - (fn + 1.0)*pi_nm1;

        const scalar w = (2.0*fn + 1.0)/(fn*(fn + 1.0));
        S1 += w*(an_[n]*pi_n  + bn_[n]*tau_n);
        S2 += w*(an_[n]*tau_n + bn_[n]*pi_n );

        // Advance pi_n to pi_{n+1}; uses recurrence index (n+1)
        // which references coefficients ((2(n+1) - 1)/((n+1) - 1) = (2n+1)/n
        // and (n+1)/n).
        if (n < N_)
        {
            const scalar pi_np1 =
                (2.0*fn + 1.0)/fn * mu * pi_n
              - (fn + 1.0)/fn * pi_nm1;
            pi_nm1 = pi_n;
            pi_n   = pi_np1;
        }
    }

    return std::norm(S1) + std::norm(S2);
}


// ************************************************************************* //
