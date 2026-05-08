/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::dragModel

Description
    Abstract base for the drag-response time provider used by the
    inertial motion model. The drag is always treated as linear in
    (U_seen - V) inside the OU exact update; the dragModel provides
    the response time τ_p* that goes into the update. Schiller-Naumann
    and other Re_p-corrected forms recover the Stokes baseline times a
    correction factor.

    Stokes baseline:
        τ_p = ρ_p d_p² / (18 μ_f)

    Re_p-corrected (Schiller-Naumann):
        τ_p* = τ_p / (1 + 0.15 Re_p^0.687),  Re_p = |U_seen - V| d_p / ν_f

    Sub-RTS family inside motionModel/inertial; not directly selectable
    from the top-level motion sub-dict.

\*---------------------------------------------------------------------------*/

#ifndef dose_dragModel_H
#define dose_dragModel_H

#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "scalar.H"
#include "vector.H"

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class dragModel Declaration
\*---------------------------------------------------------------------------*/

class dragModel
{
public:

    //- Runtime type information
    TypeName("dragModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        dragModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        dragModel() = default;


    //- Selector
    static autoPtr<dragModel> New(const dictionary& dict);


    //- Destructor
    virtual ~dragModel() = default;


    // Member Functions

        //- Drag response time τ_p* [s] for the OU update over one
        //  outer step. Inputs are physical: particle velocity V, fluid
        //  velocity seen by the particle U_seen, particle density ρ_p
        //  [kg/m³], particle diameter d_p [m], fluid kinematic
        //  viscosity ν_f [m²/s], fluid dynamic viscosity μ_f [Pa·s].
        virtual scalar tauP
        (
            const vector& V,
            const vector& Useen,
            scalar rhoP,
            scalar dP,
            scalar nuF,
            scalar muF
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
