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
    Foam::fv::opticalRadiation

Description
    fvModel wrapper for the optical-radiation library. Embed in any host
    solver via the fvModels dict to have the discrete-ordinates RTE
    solved each correction step alongside the host solver's other
    equations.

    The wrapper owns the radiance field "I" (read from the
    case time directory) and the radiationModel instance. It does NOT
    add a source term to any host equation -- the radiation field G is
    a side computation that the host can consume via the mesh registry.
    Couplings such as G driving a species growth source belong in the
    host's own fvModels.

    Example fvModels entry:
    \verbatim
        opticalRadiation
        {
            type    opticalRadiation;
            libs    ("libopticalRadiation.so");
        }
    \endverbatim

    Configuration of the radiation model itself lives in
    constant/opticalRadiationProperties as for the standalone solver.

SourceFiles
    opticalRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef fv_opticalRadiation_H
#define fv_opticalRadiation_H

#include "fvModel.H"
#include "volFieldsFwd.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class opticalRadiation Declaration
\*---------------------------------------------------------------------------*/

class opticalRadiation
:
    public fvModel
{
    // Private data

        //- Radiance placeholder field. The per-ray radiance
        //  fields are owned by the radiation model itself; this is the
        //  field used to construct the model and to satisfy the
        //  "I" entry expected at construction time.
        autoPtr<volScalarField> I_;

        //- The optical radiation model
        autoPtr<Foam::optical::radiationModel> radiationModel_;


    // Private member functions

        opticalRadiation(const opticalRadiation&) = delete;
        void operator=(const opticalRadiation&) = delete;


public:

    //- Runtime type information
    TypeName("opticalRadiation");


    // Constructors

        //- Construct from components
        opticalRadiation
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~opticalRadiation();


    // Member Functions

        // Checks

            //- This fvModel does not add source terms to host equations;
            //  it solves the RTE as a side computation and exposes G via
            //  the mesh registry.
            virtual wordList addSupFields() const;


        // Mesh changes (no-ops; the radiation model has no mesh-motion
        // or topology-change state of its own)

            virtual bool movePoints();
            virtual void topoChange(const polyTopoChangeMap&);
            virtual void mapMesh(const polyMeshMap&);
            virtual void distribute(const polyDistributionMap&);


        //- Solve the RTE. Called by the host solver each correction step.
        virtual void correct();


        //- Re-read coefficients
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
