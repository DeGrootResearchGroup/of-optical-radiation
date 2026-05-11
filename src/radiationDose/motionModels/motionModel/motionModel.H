/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::motionModel

Description
    Abstract base class for the per-particle equation of motion. The
    model owns the rule by which a particle's velocity is updated over
    one outer step, given the local fluid mean velocity U, the
    dispersion fluctuation u', and the current particle velocity.

    Two concrete models ship today:
    - tracer:   V is algebraically equal to U + u'. Recovers the
                drift-free fluid-tracer behaviour of v0.3.
    - inertial: V evolves under drag (+ optional gravity, Brownian)
                via an Ornstein-Uhlenbeck exact update. See
                inertial.H for the equations.

    The model returns two vectors per call:
    - V       : the particle's true velocity at the END of the step.
                Stored on the particle and used as V_old for the next
                step.
    - V_disp  : the mean velocity over the step (Δx/Δt). Used by the
                outer move() loop as the displacement velocity for
                trackToAndHitFace. For tracer V_disp = V; for the OU
                update V_disp = V_eq + (V_old - V_eq)·(1-exp(-ω))/ω.

    Per-track state (e.g. lagged Re_p for Schiller-Naumann) lives on
    the track object. Mirrors the dispersionModel pattern.

\*---------------------------------------------------------------------------*/

#ifndef dose_motionModel_H
#define dose_motionModel_H

#include "HashTable.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "perTrackState.H"
#include "randomGenerator.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "wordList.H"

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class motionModel Declaration
\*---------------------------------------------------------------------------*/

class motionModel
{

protected:

    // Protected Data

        const dictionary dict_;
        const fvMesh& mesh_;


public:

    // Public Classes

        //- Per-track state slot. Concrete models override with their
        //  own derived state (deriving from dose::perTrackState);
        //  models with no per-track memory use this base directly.
        //  Aliased rather than redeclared so motionModel::State and
        //  dispersionModel::State name the same polymorphic base.
        using State = perTrackState;

        //- Result of one outer-step velocity update.
        struct stepResult
        {
            vector V;       // particle velocity at end of step
            vector Vdisp;   // mean velocity over the step (for displacement)
        };


    //- Runtime type information
    TypeName("motionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        motionModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        motionModel(const dictionary& dict, const fvMesh& mesh);


    //- Selector
    static autoPtr<motionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~motionModel();


    // Static dispatch parallel to New() -- see dispersionModel.H for
    // the rationale (avoids instantiating a throwaway model just to
    // query its field requirements).

        //- Signature of a per-type requiredFields(dict) helper.
        using requiredFieldsFn = wordList(*)(const dictionary&);

        //- Look up the concrete type from `dict["type"]` and call
        //  its registered static requiredFields(dict). Fatal if no
        //  type was registered for that name.
        static wordList requiredFields(const dictionary& dict);

        //- Registration RAII. Each concrete model's .C creates one
        //  static instance to plumb itself into the dispatch table.
        struct addRequiredFields
        {
            addRequiredFields(const word& type, requiredFieldsFn fn);
        };


    // Member Functions

        //- Names of fields the model needs loaded into the registry.
        //  Default: none. Concrete models override; the override
        //  usually delegates to the model's static
        //  requiredFields(dict_).
        virtual wordList requiredFields() const
        {
            return wordList();
        }

        //- Whether the model wants a fresh unit-Gaussian sample drawn
        //  by the caller and passed in via xi. Tracer and the
        //  drag-only inertial path return false; only the Brownian-
        //  active inertial path returns true. The caller draws xi
        //  exactly once per outer step regardless of CFL re-stepping.
        virtual bool needsBrownianSample() const
        {
            return false;
        }

        //- Create a fresh per-track State object.
        virtual autoPtr<State> newState() const
        {
            return autoPtr<State>(new State());
        }

        //- Advance the particle velocity over one outer step.
        //
        //  Inputs:
        //    state   per-track state (model-owned subtype)
        //    V       particle velocity at start of step
        //    Umean   mean fluid velocity at particle position
        //    uPrime  dispersion fluctuation at particle position
        //    dt      outer-step duration
        //    xi      pre-drawn unit Gaussian (used only if
        //            needsBrownianSample() returned true)
        virtual stepResult advance
        (
            State& state,
            const vector& V,
            const vector& Umean,
            const vector& uPrime,
            scalar dt,
            const vector& xi
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
