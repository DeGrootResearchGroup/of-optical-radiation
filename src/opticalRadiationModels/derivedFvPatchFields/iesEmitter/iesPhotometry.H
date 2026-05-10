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

Class
    Foam::optical::iesPhotometry

Description
    Loader and bilinear interpolator for IESNA LM-63 Type C photometric
    files. The candela table is held in memory and queried on demand by
    the iesEmitter boundary condition.

    Only Type C photometry is supported. Type C convention: vertical
    angle gamma is measured from the fixture's nadir (gamma=0 along the
    nominal beam axis); horizontal angle h is measured around that
    axis. Both are in DEGREES on disk and at the interpolate() API.

    Symmetry of the horizontal angle list is inferred from the table:
      - 1 angle -> ROTATIONAL (axisymmetric about the beam axis)
      - last == 90 -> QUADRANT (4-fold mirror symmetry)
      - last == 180 -> BILATERAL (mirror about h=180)
      - otherwise (typically last == 360) -> FULL (no implicit symmetry)

    Bilinear interpolation in (gamma, h) is used; out-of-range
    horizontal queries are folded back into the stored range using the
    inferred symmetry. Out-of-range vertical queries are clamped at the
    table edges.

    Values in the table are returned scaled by the IES candela
    multiplier; the iesEmitter BC further renormalises against a
    user-supplied total radiant flux, so absolute photometric units
    (candela vs W/sr) do not matter.

SourceFiles
    iesPhotometry.C

\*---------------------------------------------------------------------------*/

#ifndef optical_iesPhotometry_H
#define optical_iesPhotometry_H

#include "scalarList.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                       Class iesPhotometry Declaration
\*---------------------------------------------------------------------------*/

class iesPhotometry
{
public:

    enum hSymmetry
    {
        ROTATIONAL,
        QUADRANT,
        BILATERAL,
        FULL
    };

private:

    //- Source file (kept for diagnostics + write())
    fileName path_;

    //- Vertical angles [degrees], strictly ascending
    scalarList vertical_;

    //- Horizontal angles [degrees], strictly ascending
    scalarList horizontal_;

    //- Candela values, indexed candela_[ih*nVertical + iv]
    //  (intensity at horizontal_[ih], vertical_[iv]); already
    //  multiplied by the IES candela multiplier.
    scalarList candela_;

    //- Inferred horizontal symmetry
    hSymmetry symmetry_;


    // Private member functions

        //- Parse the file at construction
        void load_();

        //- Infer hSymmetry_ from horizontal angle list
        void inferSymmetry_();

        //- Map a horizontal angle in [0, 360) into the stored range using
        //  the inferred symmetry
        scalar foldHorizontal_(scalar h) const;

public:

    // Constructors

        //- Construct from a file path
        explicit iesPhotometry(const fileName& path);


    // Member functions

        //- Bilinear interpolation. Inputs in DEGREES.
        //- gamma is clamped into [vertical_.first(), vertical_.last()].
        //- h is wrapped to [0, 360) and folded by symmetry.
        scalar interpolate(scalar gammaDeg, scalar hDeg) const;

        //- Source path
        const fileName& path() const { return path_; }

        //- Number of vertical angles
        label nVertical() const { return vertical_.size(); }

        //- Number of horizontal angles in the stored table (pre-symmetry)
        label nHorizontal() const { return horizontal_.size(); }

        //- Inferred horizontal symmetry
        hSymmetry symmetry() const { return symmetry_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
