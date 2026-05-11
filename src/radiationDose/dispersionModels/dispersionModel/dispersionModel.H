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
#include "dictionary.H"
#include "fvMesh.H"
#include "perTrackState.H"
#include "randomGenerator.H"
#include "runTimeSelectionTables.H"
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

        //- Velocity fluctuation u' [m/s] at position x in cell celli
        //  over a step of size dt [s]. The State reference is the
        //  track's own per-particle state (not shared across tracks),
        //  so concurrent calls from different threads are race-free.
        //  Concrete models look up their own turbulence fields
        //  (k, epsilon, ...) from the registry as needed.
        virtual vector fluctuation
        (
            State& state,
            const vector& x,
            label celli,
            scalar dt,
            randomGenerator& rng
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
