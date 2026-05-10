/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

#include "iesPhotometry.H"
#include "error.H"

#include <fstream>
#include <string>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::iesPhotometry::iesPhotometry(const fileName& path)
:
    path_(path),
    vertical_(0),
    horizontal_(0),
    candela_(0),
    symmetry_(ROTATIONAL)
{
    load_();
    inferSymmetry_();
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::optical::iesPhotometry::load_()
{
    std::ifstream is(path_.c_str());
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot open IES file: " << path_
            << exit(FatalError);
    }

    // Phase 1: line-based scan for the TILT= marker. Everything before
    // it -- version line, [KEYWORD] metadata, blank lines -- is
    // ignored. The IES spec requires a TILT line; if it's missing we
    // can't trust the file.
    std::string line;
    std::string tiltMode;
    bool foundTilt = false;
    while (std::getline(is, line))
    {
        const auto pos = line.find("TILT=");
        if (pos != std::string::npos)
        {
            tiltMode = line.substr(pos + 5);
            const auto first = tiltMode.find_first_not_of(" \t\r\n");
            if (first == std::string::npos)
            {
                tiltMode.clear();
            }
            else
            {
                const auto last = tiltMode.find_last_not_of(" \t\r\n");
                tiltMode = tiltMode.substr(first, last - first + 1);
            }
            foundTilt = true;
            break;
        }
    }
    if (!foundTilt)
    {
        FatalErrorInFunction
            << "TILT= line not found in IES file: " << path_
            << exit(FatalError);
    }

    // Phase 2: token-based reads. From here on the format is purely
    // whitespace-separated numerical fields, possibly spanning many
    // lines. operator>> handles both.
    if (tiltMode != "NONE")
    {
        // KrCl downlights typically use TILT=NONE. Fail loud rather
        // than silently miscompute the angular distribution.
        FatalErrorInFunction
            << "Only TILT=NONE is supported; got TILT=" << tiltMode
            << " in IES file: " << path_
            << exit(FatalError);
    }

    // The 10 fixed header values
    label numLamps;
    scalar lumensPerLamp, candelaMultiplier;
    label nV, nH, photometricType, unitsType;
    scalar lumOpenW, lumOpenL, lumOpenH;
    is >> numLamps
       >> lumensPerLamp
       >> candelaMultiplier
       >> nV
       >> nH
       >> photometricType
       >> unitsType
       >> lumOpenW
       >> lumOpenL
       >> lumOpenH;

    if (!is.good())
    {
        FatalErrorInFunction
            << "Failed reading header block from IES file: " << path_
            << exit(FatalError);
    }

    if (photometricType != 1)
    {
        FatalErrorInFunction
            << "Only Type C photometry (photometricType = 1) is"
            << " supported; got photometricType = " << photometricType
            << " in IES file: " << path_
            << exit(FatalError);
    }

    if (nV < 2)
    {
        FatalErrorInFunction
            << "IES file must have at least 2 vertical angles; got "
            << nV << " in " << path_
            << exit(FatalError);
    }

    if (nH < 1)
    {
        FatalErrorInFunction
            << "IES file must have at least 1 horizontal angle; got "
            << nH << " in " << path_
            << exit(FatalError);
    }

    // Ballast block (3 values: ballastFactor, future use, input watts)
    scalar ballastFactor, balLampPhotoFactor, inputWatts;
    is >> ballastFactor >> balLampPhotoFactor >> inputWatts;

    // Vertical angles
    vertical_.setSize(nV);
    for (label i = 0; i < nV; ++i)
    {
        is >> vertical_[i];
    }
    for (label i = 1; i < nV; ++i)
    {
        if (vertical_[i] <= vertical_[i - 1])
        {
            FatalErrorInFunction
                << "Vertical angles must be strictly ascending in "
                << path_ << "; got vertical_[" << i - 1 << "] = "
                << vertical_[i - 1] << ", vertical_[" << i << "] = "
                << vertical_[i] << exit(FatalError);
        }
    }

    // Horizontal angles
    horizontal_.setSize(nH);
    for (label i = 0; i < nH; ++i)
    {
        is >> horizontal_[i];
    }
    for (label i = 1; i < nH; ++i)
    {
        if (horizontal_[i] <= horizontal_[i - 1])
        {
            FatalErrorInFunction
                << "Horizontal angles must be strictly ascending in "
                << path_ << "; got horizontal_[" << i - 1 << "] = "
                << horizontal_[i - 1] << ", horizontal_[" << i << "] = "
                << horizontal_[i] << exit(FatalError);
        }
    }

    // Candela values: nH blocks of nV values each.
    // Storage layout matches: candela_[ih*nV + iv].
    //
    // The post-loop `is.good() && !is.eof()` check we used to do was
    // effectively a no-op: `>>` on a truncated stream sets failbit AND
    // eofbit, so the conjunction is always false. Instead check
    // is.fail() after each read so a truncated or malformed file
    // fatals at the exact value where parsing stopped, rather than
    // silently leaving the remainder of `candela_` at uninitialised /
    // zero values.
    candela_.setSize(nH*nV);
    for (label k = 0; k < nH*nV; ++k)
    {
        is >> candela_[k];
        if (is.fail())
        {
            const label ih = k/nV;
            const label iv = k - ih*nV;
            FatalErrorInFunction
                << "Failed reading candela value at horizontal index "
                << ih << ", vertical index " << iv
                << " (flat position " << k << " of " << nH*nV
                << ") from IES file: " << path_
                << "; file may be truncated or contain non-numeric data"
                << exit(FatalError);
        }
    }

    // Apply candela multiplier and ballast factor (the latter is
    // typically 1 for non-fluorescent fixtures).
    const scalar scale = candelaMultiplier*ballastFactor;
    forAll(candela_, i)
    {
        candela_[i] *= scale;
        if (candela_[i] < 0.0)
        {
            // IES values are non-negative by definition; clip
            // numerical negatives at zero rather than propagating
            // them into the BC's renormalisation.
            candela_[i] = 0.0;
        }
    }
}


void Foam::optical::iesPhotometry::inferSymmetry_()
{
    const label nH = horizontal_.size();
    if (nH == 1)
    {
        symmetry_ = ROTATIONAL;
        return;
    }

    const scalar tol = 1e-3;
    const scalar last = horizontal_.last();

    // Standard IES Type C horizontal-range conventions:
    //   first == 0 always (we don't enforce it but it's implicit).
    //   last  ~= 90  -> 4-fold symmetry (quadrant)
    //   last  ~= 180 -> bilateral symmetry (mirror about h=180)
    //   last  ~= 360 (or matches first) -> full table.
    if (mag(last - 90.0) < tol)
    {
        symmetry_ = QUADRANT;
    }
    else if (mag(last - 180.0) < tol)
    {
        symmetry_ = BILATERAL;
    }
    else
    {
        symmetry_ = FULL;
    }
}


Foam::scalar Foam::optical::iesPhotometry::foldHorizontal_(scalar h) const
{
    // Caller ensures h in [0, 360).
    switch (symmetry_)
    {
        case ROTATIONAL:
            return horizontal_[0];

        case QUADRANT:
            // 4-fold: [0..90] is the master quadrant.
            if (h > 270.0) return 360.0 - h;
            if (h > 180.0) return h - 180.0;
            if (h >  90.0) return 180.0 - h;
            return h;

        case BILATERAL:
            // 2-fold: [0..180] master.
            if (h > 180.0) return 360.0 - h;
            return h;

        case FULL:
        default:
            return h;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::optical::iesPhotometry::interpolate
(
    scalar gammaDeg,
    scalar hDeg
) const
{
    const label nV = vertical_.size();
    const label nH = horizontal_.size();

    // Wrap horizontal into [0, 360), then fold by symmetry into the
    // stored range.
    while (hDeg >= 360.0) hDeg -= 360.0;
    while (hDeg <    0.0) hDeg += 360.0;
    const scalar hFolded = foldHorizontal_(hDeg);

    // Vertical bracket. Clamp into the stored range -- IES tables are
    // not extrapolated beyond the published gamma range; the boundary
    // values are the best available estimate.
    label iv;
    scalar fv;
    if (gammaDeg <= vertical_.first())
    {
        iv = 0;
        fv = 0.0;
    }
    else if (gammaDeg >= vertical_.last())
    {
        iv = nV - 2;
        fv = 1.0;
    }
    else
    {
        iv = 0;
        while (iv < nV - 2 && vertical_[iv + 1] < gammaDeg) ++iv;
        fv =
            (gammaDeg - vertical_[iv])
           /(vertical_[iv + 1] - vertical_[iv]);
    }

    // 1-D interpolate (vertical only) for axisymmetric tables.
    if (nH == 1)
    {
        return
            (1.0 - fv)*candela_[0*nV + iv]
           +       fv *candela_[0*nV + iv + 1];
    }

    // FULL-symmetry seam wrap. For ROTATIONAL/QUADRANT/BILATERAL,
    // foldHorizontal_ folds hDeg into [horizontal_.first(),
    // horizontal_.last()] and the standard bracket below works. For
    // FULL, foldHorizontal_ is a no-op; if the table closes early
    // (e.g. last horizontal angle at 350 deg instead of 360 deg) then
    // hFolded can fall in the cyclic [last, 360 + first) gap, and the
    // standard bracket would extrapolate beyond the last column.
    // Handle the gap explicitly by interpolating cyclically between
    // candela_[ih = nH - 1] and candela_[ih = 0], with a horizontal
    // span of (first + 360 - last).
    if (symmetry_ == FULL && hFolded > horizontal_.last())
    {
        const scalar h0 = horizontal_.last();
        const scalar h1 = horizontal_.first() + 360.0;
        const scalar fh = (hFolded - h0)/(h1 - h0);
        const label ih0 = nH - 1;
        const label ih1 = 0;
        const scalar I00 = candela_[ih0*nV + iv    ];
        const scalar I01 = candela_[ih0*nV + iv + 1];
        const scalar I10 = candela_[ih1*nV + iv    ];
        const scalar I11 = candela_[ih1*nV + iv + 1];
        return
            (1.0 - fh)*((1.0 - fv)*I00 + fv*I01)
          +        fh *((1.0 - fv)*I10 + fv*I11);
    }

    // Horizontal bracket. foldHorizontal_ guarantees hFolded lies in
    // [horizontal_.first(), horizontal_.last()] for the
    // ROTATIONAL/QUADRANT/BILATERAL cases; the FULL seam-wrap case
    // above handles the hFolded > horizontal_.last() shortfall.
    label ih = 0;
    while (ih < nH - 2 && horizontal_[ih + 1] < hFolded) ++ih;
    const scalar fh =
        (hFolded - horizontal_[ih])
       /(horizontal_[ih + 1] - horizontal_[ih]);

    const scalar I00 = candela_[ih      *nV + iv    ];
    const scalar I01 = candela_[ih      *nV + iv + 1];
    const scalar I10 = candela_[(ih + 1)*nV + iv    ];
    const scalar I11 = candela_[(ih + 1)*nV + iv + 1];

    return
        (1.0 - fh)*((1.0 - fv)*I00 + fv*I01)
      +        fh *((1.0 - fv)*I10 + fv*I11);
}


// ************************************************************************* //
