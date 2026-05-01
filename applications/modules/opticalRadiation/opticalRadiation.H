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
    Foam::solvers::opticalRadiation

Description
    Solver module for the optical-radiation library. Designed to be loaded
    per region by `foamMultiRun` so that a multi-region case can solve
    radiation in each region while a `refractiveCoupled`
    boundary condition handles the cross-region coupling at refractive
    interfaces.

    For single-region cases the standalone `opticalRadiationFoam` binary
    is simpler. For regions that need radiation alongside another
    primary physics (flow, solid heat conduction), use the
    `Foam::fv::opticalRadiation` fvModel wrapper instead, embedded
    inside the host's solver module.

    Example controlDict for a two-region case driven by `foamMultiRun`:
    \verbatim
        regionSolvers
        {
            mediumA   opticalRadiation;
            mediumB   opticalRadiation;
        }
    \endverbatim

SourceFiles
    opticalRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef solvers_opticalRadiation_H
#define solvers_opticalRadiation_H

#include "solver.H"
#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                       Class opticalRadiation Declaration
\*---------------------------------------------------------------------------*/

class opticalRadiation
:
    public solver
{
    // Region-local owned data

        //- Radiation intensity field "I" (read from this region's 0/I).
        autoPtr<volScalarField> I_;

        //- The optical radiation model for this region.
        autoPtr<Foam::optical::radiationModel> radiationModel_;


    // Disallow copy / assignment
    opticalRadiation(const opticalRadiation&) = delete;
    void operator=(const opticalRadiation&) = delete;


public:

    //- Runtime type information
    TypeName("opticalRadiation");


    // Constructors

        //- Construct from region mesh
        opticalRadiation(fvMesh& mesh);


    //- Destructor
    virtual ~opticalRadiation();


    // Member Functions

        //- No CFL-style time-step constraint from radiation
        virtual scalar maxDeltaT() const;

        //- Solve the RTE once at the start of the time-step.
        //  preSolve() is called before moveMesh(), so this is correct for
        //  static meshes only. For moving-mesh cases the radiation solve
        //  should be moved to prePredictor() (called once per outer PIMPLE
        //  iteration, after motion); this is not currently used so kept
        //  here for simplicity.
        virtual void preSolve();


        // --- No-ops: pure-virtuals from `solver` that don't apply to a
        //     pure-radiation solver module. (The host PIMPLE loop calls
        //     these on every region; for radiation they are unused.)

        virtual void moveMesh()                             {}
        virtual void motionCorrector()                      {}
        virtual void prePredictor()                         {}
        virtual void momentumTransportPredictor()           {}
        virtual void thermophysicalTransportPredictor()     {}
        virtual void momentumPredictor()                    {}
        virtual void thermophysicalPredictor()              {}
        virtual void pressureCorrector()                    {}
        virtual void momentumTransportCorrector()           {}
        virtual void thermophysicalTransportCorrector()     {}
        virtual void postSolve()                            {}


        //- regIOobject I/O. Per-region radiation state (I, GLambda_*, etc.)
        //  is written via each field's own AUTO_WRITE flag, so the solver
        //  itself has nothing extra to serialise here.
        virtual bool writeData(Ostream& os) const;

        //- Re-read radiation properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
