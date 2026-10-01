/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  | Copyright (C) 2018-2026 DeGroot Research Group
-------------------------------------------------------------------------------

Class
    Foam::dose::dispersionModel

Description
    Abstract base class for stochastic turbulent dispersion models.
    Returns a velocity fluctuation u' to add to the mean velocity at
    a given particle position.

    One outer step talks to the model in three calls, because the
    step's duration is not known until the fluctuation is: the CFL
    bound depends on U + u'.

    - beginStep(): sample the local turbulence at the step's start and
      draw every random number the step may need, given the longest
      duration it can have (dtMax).
    - fluctuation(): the mean u' over a step of a given duration. Pure;
      the tracker calls it once for dtMax and again for the CFL-reduced
      duration, and both see the same random draws. The CFL bound
      itself uses cflFluctuation() and maxStep() (see below).
    - endStep(): advance the per-track state by the duration actually
      taken.

    correct() rebuilds whatever the model caches from the registry's
    turbulence fields; the function object calls it before each run,
    so a host solver that updates the fields between calls is seen.

    Per-track state (e.g. DRW eddy memory) lives on the track object
    rather than inside the model. The model defines a State subclass
    plus a newState() factory; the track owns one State pointer.
    Eliminating the model's shared state map makes the integration
    loop trivially thread-parallel.

\*---------------------------------------------------------------------------*/

#ifndef dose_dispersionModel_H
#define dose_dispersionModel_H

#include "HashTable.H"
#include "autoPtr.H"
#include "barycentric.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "perTrackState.H"
#include "randomGenerator.H"
#include "runTimeSelectionTables.H"
#include "tetIndices.H"
#include "vector.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class dispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersionModel
{

protected:

    // Protected Data

        //- Model dictionary (the contents of the "dispersion" sub-dict)
        const dictionary dict_;

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    // Public Classes

        //- Per-track state slot. Each concrete model defines its own
        //  derived state (deriving from dose::perTrackState, e.g.
        //  discreteRandomWalk::DRWState carries the current eddy
        //  fluctuation and remaining lifetime). Models with no per-
        //  track state use this base directly. Aliased rather than
        //  redeclared so motionModel::State and dispersionModel::State
        //  name the same polymorphic base.
        using State = perTrackState;


    //- Runtime type information
    TypeName("dispersionModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dispersionModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        dispersionModel(const dictionary& dict, const fvMesh& mesh);


    //- Selector
    static autoPtr<dispersionModel> New
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~dispersionModel();


    // Static dispatch parallel to New() -- inspecting the model's
    // field requirements without instantiating it. The instance form
    // of requiredFields() below is still available for callers that
    // already have a model; the static form is the cheaper path for
    // a function-object that only needs to know which fields to
    // request from the host before building anything.

        //- Signature of a per-type requiredFields(dict) helper.
        using requiredFieldsFn = wordList(*)(const dictionary&);

        //- Look up the concrete type from `dict["type"]` and call
        //  its registered static requiredFields(dict). Fatal if no
        //  type was registered for that name.
        static wordList requiredFields(const dictionary& dict);

        //- Registration RAII. Each concrete model's .C creates one
        //  static instance to plumb itself into requiredFieldsTable_.
        struct addRequiredFields
        {
            addRequiredFields(const word& type, requiredFieldsFn fn);
        };


    // Member Functions

        //- Names of fields the model needs loaded into the registry
        //  (e.g. k, epsilon for DRW). Default: none. Concrete models
        //  override; the override usually delegates to the model's
        //  static requiredFields(dict_).
        virtual wordList requiredFields() const
        {
            return wordList();
        }

        //- Create a fresh per-track State object. The track owns it.
        //  Default: a base State with no per-track memory.
        virtual autoPtr<State> newState() const
        {
            return autoPtr<State>(new State());
        }

        //- Rebuild anything cached from the registry's fields (e.g.
        //  interpolators). Called by the function object before each
        //  run, after the fields it reads are current. Default: none.
        virtual void correct()
        {}

        //- Start an outer step at the particle's position: sample the
        //  turbulence there and draw the random numbers a step of up
        //  to dtMax [s] needs. Default: nothing to do.
        virtual void beginStep
        (
            State& state,
            const barycentric& coordinates,
            const tetIndices& tetIs,
            scalar dtMax,
            randomGenerator& rng
        ) const
        {}

        //- Mean velocity fluctuation u' [m/s] over a step of duration
        //  dt <= dtMax [s] begun by beginStep(). Reads the state only,
        //  so it may be called for several candidate durations. The
        //  State is the track's own (never shared across tracks), so
        //  concurrent calls from different threads are race-free.
        virtual vector fluctuation(const State& state, scalar dt) const = 0;

        //- The fluctuation the tracker's CFL bound sees for a step of
        //  duration dt. Default: fluctuation() itself. A model whose
        //  step displacement is a fresh random draw (randomDisplacement)
        //  returns its deterministic part instead: bounding the step
        //  by the draw it takes would make large draws take short
        //  steps, and the walk would spread too little.
        virtual vector cflFluctuation(const State& state, scalar dt) const
        {
            return fluctuation(state, dt);
        }

        //- The longest step [s] the model allows, given the largest
        //  displacement maxDisplacement [m] a step may make (cflMax
        //  times the cell size). Default: no limit.
        virtual scalar maxStep
        (
            const State& state,
            scalar maxDisplacement
        ) const
        {
            return great;
        }

        //- Finish the step: advance the state by the duration dt [s]
        //  actually taken. Default: nothing to do.
        virtual void endStep(State& state, scalar dt) const
        {}

        //- The particle reflected off a boundary face with unit normal
        //  n: reflect whatever velocity the state carries into the
        //  rest of the step and the next one. Default: nothing to do.
        virtual void reflect(State& state, const vector& n) const
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
