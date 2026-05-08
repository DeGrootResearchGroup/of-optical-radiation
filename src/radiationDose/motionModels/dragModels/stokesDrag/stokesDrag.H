/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::stokesDrag

Description
    Stokes drag (low-Re limit). Returns the constant response time

        τ_p = ρ_p d_p² / (18 μ_f)

    irrespective of |U_seen - V|. Valid strictly for Re_p ≪ 1; for
    larger Re_p select schillerNaumann instead.

\*---------------------------------------------------------------------------*/

#ifndef dose_stokesDrag_H
#define dose_stokesDrag_H

#include "dragModel.H"

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class stokesDrag Declaration
\*---------------------------------------------------------------------------*/

class stokesDrag
:
    public dragModel
{
public:

    //- Runtime type information
    TypeName("stokes");


    // Constructors

        stokesDrag(const dictionary& dict);


    //- Destructor
    virtual ~stokesDrag() = default;


    // Member Functions

        virtual scalar tauP
        (
            const vector& /*V*/,
            const vector& /*Useen*/,
            scalar rhoP,
            scalar dP,
            scalar /*nuF*/,
            scalar muF
        ) const
        {
            return rhoP*dP*dP/(18.0*muF);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
