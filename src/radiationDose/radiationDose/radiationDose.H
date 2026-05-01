/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | radiationDose: Lagrangian radiation dose tracking
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::functionObjects::radiationDose

Description
    Function object that integrates radiation dose along Lagrangian
    particle trajectories through a (typically frozen) flow field.

    Reads a velocity field U [m/s] and a fluence-rate field G [W/m^2],
    seeds particles through a configurable seeding model, integrates
    their motion (with optional stochastic turbulent dispersion), and
    accumulates dose D along each trajectory.

    Units convention:
    - Fluence rate input G is in W/m^2 (OpenFOAM SI standard).
    - Accumulated dose D is in mJ/cm^2 (UV reactor literature standard).
    - The unit conversion factor 0.1 (W/m^2 . s -> mJ/cm^2) is applied
      inside the integration kernel.
    - Inactivation rate kInact is in cm^2/mJ.

SourceFiles
    radiationDose.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_radiationDose_H
#define functionObjects_radiationDose_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class radiationDose Declaration
\*---------------------------------------------------------------------------*/

class radiationDose
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of the velocity field [m/s]
        word UName_;

        //- Name of the fluence-rate field [W/m^2]
        word GName_;


public:

    //- Conversion factor: G [W/m^2] * dt [s] -> D [mJ/cm^2]
    //  1 J/m^2 = 0.1 mJ/cm^2 (since 1 m^2 = 1e4 cm^2 and 1 J = 1e3 mJ)
    static constexpr scalar Wm2_s_to_mJcm2 = 0.1;


    //- Runtime type information
    TypeName("radiationDose");


    // Constructors

        //- Construct from Time and dictionary
        radiationDose
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        radiationDose(const radiationDose&) = delete;


    //- Destructor
    virtual ~radiationDose();


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute (integrate trajectories and accumulate dose)
        virtual bool execute();

        //- Write outputs (VTK tracks, dose CSV, summary)
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const radiationDose&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
