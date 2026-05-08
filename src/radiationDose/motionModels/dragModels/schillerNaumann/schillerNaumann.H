/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::schillerNaumann

Description
    Schiller-Naumann drag correction:

        τ_p* = τ_p / (1 + 0.15 Re_p^0.687)
        Re_p = |U_seen - V| d_p / ν_f

    Reduces to Stokes for Re_p → 0, valid up to Re_p ≈ 1000. Re_p is
    evaluated from the velocities passed in (no lagging — the OU
    update treats the result as the linear-drag coefficient over one
    outer step, so a single Re_p evaluation per step is consistent
    with the linearisation).

\*---------------------------------------------------------------------------*/

#ifndef dose_schillerNaumann_H
#define dose_schillerNaumann_H

#include "dragModel.H"

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class schillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class schillerNaumann
:
    public dragModel
{
public:

    //- Runtime type information
    TypeName("schillerNaumann");


    // Constructors

        schillerNaumann(const dictionary& dict);


    //- Destructor
    virtual ~schillerNaumann() = default;


    // Member Functions

        virtual scalar tauP
        (
            const vector& V,
            const vector& Useen,
            scalar rhoP,
            scalar dP,
            scalar nuF,
            scalar muF
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
