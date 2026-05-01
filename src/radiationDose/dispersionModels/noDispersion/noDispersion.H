/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::dose::noDispersion

Description
    Trivial dispersion model that returns zero velocity fluctuation
    (deterministic streamlines).

\*---------------------------------------------------------------------------*/

#ifndef dose_noDispersion_H
#define dose_noDispersion_H

#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dose
{

/*---------------------------------------------------------------------------*\
                       Class noDispersion Declaration
\*---------------------------------------------------------------------------*/

class noDispersion
:
    public dispersionModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        noDispersion(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~noDispersion() = default;


    // Member Functions

        virtual vector fluctuation
        (
            label trackId,
            const vector& x,
            label celli,
            scalar dt,
            randomGenerator& rng
        ) const
        {
            return vector::zero;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dose
} // End namespace Foam

#endif

// ************************************************************************* //
