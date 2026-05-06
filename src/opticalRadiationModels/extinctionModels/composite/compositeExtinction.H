/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::optical::compositeExtinction

Description
    Sums the absorption and scattering coefficients of an arbitrary set of
    child extinction models. Child models are listed by user-chosen name
    inside compositeCoeffs.models; each entry is the same dictionary an
    extinction model would consume at the top level.

    Example: combine Rayleigh scattering with a species-dependent absorber:

        extinctionModel composite;

        compositeCoeffs
        {
            nBands 1;

            models
            {
                rayleigh
                {
                    extinctionModel rayleigh;
                    rayleighCoeffs
                    {
                        nBands       1;
                        wavelengths  (222);
                        T            293.15;
                        p            101325;
                    }
                }

                ozoneAbsorption
                {
                    extinctionModel linearSpecies;
                    linearSpeciesCoeffs
                    {
                        nBands           1;
                        absorption       true;
                        nAbsorbing       1;
                        absorbingVars    (O3);
                        absorptionCoeff  ((1310));   // 222 nm O3, m^2/mol
                        scattering       false;
                        nScattering      0;
                    }
                }
            }
        }

SourceFiles
    compositeExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef opticalRadiationCompositeExtinction_H
#define opticalRadiationCompositeExtinction_H

#include "extinctionModel.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                  Class compositeExtinction Declaration
\*---------------------------------------------------------------------------*/

class compositeExtinction
:
    public extinctionModel
{

    // Private data

        //- Extinction model dictionary
        dictionary coeffsDict_;

        //- Number of bands
        label nBands_;

        //- Child extinction models
        PtrList<extinctionModel> models_;

public:

    //- Runtime type information
    TypeName("composite");


    // Constructors

        //- Construct from components
        compositeExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~compositeExtinction();


    // Member Functions

        //- Update children, then sum into ALambda_/SLambda_
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
