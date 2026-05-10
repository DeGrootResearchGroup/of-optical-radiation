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

Class
    Foam::optical::phaseFunctionModel

Description
    Base class for opticalRadiation scattering phase functions.

    The discrete in-scatter source assembled in DOM::calculate is

        ds(rayI) = sum over j of (I_j * table[i, j, iBand])

    where the table entries are bin integrals of the angular phase
    function Phi(s_i . s'):

        table[i, j, iBand] ~ integral over Omega_j of
                             Phi(s_i . s') dOmega'

    (suitably normalised; see buildPhaseTable below).

    Concrete subclasses implement only the angular shape phaseShape(cosV,
    iBand); buildPhaseTable handles the Murthy-Mathur sub-pixel sampling
    over each Omega_j bin and the row-wise normalisation.

SourceFiles
    phaseFunctionModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseFunctionModel_H
#define phaseFunctionModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{
class DOM;

/*---------------------------------------------------------------------------*\
                       Class phaseFunctionModel Declaration
\*---------------------------------------------------------------------------*/

class phaseFunctionModel
{
protected:

    // Protected data

        //- Owning DOM (for ray angles, dimensions, etc.)
        const DOM& dom_;

        //- opticalRadiation model dictionary
        const dictionary dict_;

        //- Number of spatial dimensions (2 or 3); selects pixel layout in
        //  buildPhaseTable.
        label nDim_;

        //- Number of spectral bands.
        label nBand_;

        //- Number of angular bins per band (= 2*nPhi*nTheta).
        label nAngle_;

        //- Set true when the concrete model contributes to the in-scatter
        //  source. Subclass constructor reads this from its own
        //  ${modelType}Coeffs sub-dict.
        bool inScatter_;

        //- Number of sub-divisions per pixel side in the table-construction
        //  pixel sum (subAngleNum^2 sub-pixels per Omega_j bin in 3-D,
        //  subAngleNum in 2-D).
        label subAngleNum_;

        //- Flat-indexed phase function table; entry [i, j, iBand] lives at
        //  phaseFunction_[j + i*nAngle_ + iBand*nAngle_*nAngle_]. Filled by
        //  buildPhaseTable; consumed by correct().
        DynamicList<scalar, 1, 0> phaseFunction_;


    // Protected member functions

        //- Angular shape Phi(cosV) of the phase function. Subclasses
        //  override; absolute scaling does not matter (the row-
        //  normalisation in buildPhaseTable absorbs it). The default
        //  returns 1.0 so a subclass that hasn't overridden produces an
        //  isotropic table -- this is also what the base class itself
        //  evaluates to when instantiated directly via the missing-
        //  dictionary fallback in phaseFunctionModel::New() (in which
        //  case buildPhaseTable is not called at all and the in-scatter
        //  source is skipped because inScatter_ defaults to false).
        virtual scalar phaseShape
        (
            const scalar cosV,
            const label iBand
        ) const
        {
            return 1.0;
        }

        //- Fill phaseFunction_ from phaseShape(). Subclasses call this
        //  from their constructor *after* reading any per-band parameters
        //  (so virtual dispatch reaches the derived phaseShape).
        //
        //  Algorithm per (iBand, source ray i):
        //
        //    1. For each destination ray j, sub-sample Omega_j on a regular
        //       (subAngleNum)-by-(subAngleNum) pixel grid (3-D) or
        //       subAngleNum-strip (2-D), and accumulate the bin integral
        //       Psi_ij = sum over pixels of phaseShape(s_i . s_pixel) *
        //                                  d Omega_pixel
        //       which approximates integral over Omega_j of Phi(s_i . s')
        //       dOmega'. The pixel weight d Omega_pixel uses the standard
        //       sin(theta) Jacobian in 3-D and 2*d phi in 2-D.
        //
        //    2. Row-normalise: divide every Psi_ij in row i by their sum
        //       pfSum_i = sum over j of Psi_ij. The result satisfies
        //       sum over j of table[i, j, iBand] = 1. Without this step,
        //       the absolute scaling of phaseShape leaks into the
        //       in-scatter source magnitude (problematic for shapes that
        //       are not normalised to integrate to 1 over 4 pi, e.g. raw
        //       BHMIE |S1|^2 + |S2|^2 from mieKernel::phaseIntensity).
        //
        //  The runtime in-scatter sum in DOM::calculate then computes
        //  ds(rayI) = sum over j of (I_j * table[i, j, iBand]) directly,
        //  with no extra solid-angle weight.
        void buildPhaseTable();


private:

    // Private member functions

        phaseFunctionModel(const phaseFunctionModel&) = delete;
        void operator=(const phaseFunctionModel&) = delete;


public:

    //- Runtime type information
    TypeName("phaseFunctionModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseFunctionModel,
        dictionary,
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        ),
        (dom, dict, nDim)
    );


    // Constructors

        //- Construct from components.
        phaseFunctionModel
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    // Selector

        static autoPtr<phaseFunctionModel> New
        (
            const DOM& dom,
            const dictionary& dict,
            const label& nDim
        );


    //- Destructor
    virtual ~phaseFunctionModel();


    // Member functions

        // Access

            //- Reference to the dictionary the model was constructed from.
            inline const dictionary& dict() const
            {
                return dict_;
            }

            //- Whether this model contributes to the in-scatter source.
            virtual bool inScatter() const
            {
                return inScatter_;
            }

            //- Phase-function table value table[rayI, rayJ, iBand]. The
            //  rayI / rayJ flat indices already include the band offset
            //  (rayI = iAngleI + iBand*nAngle_); the table is laid out by
            //  per-band angle indices, so we strip the offset before
            //  forming the lookup.
            virtual scalar correct
            (
                const label rayI,
                const label rayJ,
                const label iBand
            ) const;

            //- Pointer to the contiguous row table[rayI, *, iBand] of
            //  length nAngle_. Lets the in-scatter loop in DOM hoist
            //  the (rayI, iBand) lookup out of the inner ray-pair loop
            //  -- one virtual call per (rayI, iBand) pair instead of
            //  one per (rayI, rayJ, iBand) tuple.
            //
            //  Returns nullptr if the table has not been built (e.g.
            //  inScatter == false, or this is the no-op base class
            //  used as the missing-dictionary fallback).
            const scalar* phaseRow
            (
                const label rayI,
                const label iBand
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
