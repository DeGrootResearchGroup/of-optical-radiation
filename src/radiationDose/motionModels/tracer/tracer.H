/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::tracer

Description
    Trivial fluid-tracer motion model. The particle velocity is
    algebraically equal to the local fluid mean plus the dispersion
    fluctuation:

        V = U + u'

    Both the end-of-step velocity and the displacement velocity equal
    this sum (no integration). Recovers the v0.3 fluid-tracer behaviour
    when selected as `motion { type tracer; }` (or omitted, which
    defaults to tracer).

\*---------------------------------------------------------------------------*/

#ifndef dose_tracer_H
#define dose_tracer_H

#include "motionModel.H"

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class tracer Declaration
\*---------------------------------------------------------------------------*/

class tracer
:
    public motionModel
{
public:

    //- Runtime type information
    TypeName("tracer");


    // Constructors

        tracer(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~tracer() = default;


    // Member Functions

        virtual stepResult advance
        (
            State& /*state*/,
            const vector& /*V*/,
            const vector& Umean,
            const vector& uPrime,
            scalar /*dt*/,
            const vector& /*xi*/
        ) const
        {
            const vector Vseen = Umean + uPrime;
            return {Vseen, Vseen};
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
