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
    Foam::optical::linearSpeciesExtinction

Description
    Model for constant absorption and scattering coefficients

SourceFiles
    linearSpeciesExtinction.C

\*---------------------------------------------------------------------------*/

#ifndef opticalRadiationWideBandExtinction_H
#define opticalRadiationWideBandExtinction_H

#include "extinctionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace optical
{

/*---------------------------------------------------------------------------*\
                  Class linearSpeciesExtinction Declaration
\*---------------------------------------------------------------------------*/

class linearSpeciesExtinction
:
    public extinctionModel
{
	
    // Private data
        
        //- Extinction model dictionary
        dictionary coeffsDict_;
        
        //- Number of bands
        label nBands_;

        //- Switch for absorption model
        Switch absorption_;

        //- Number of absorbing species
        label nAbsorbing_;

        //- Names of absorbing species
        DynamicList<word, 1, 0> aSpecies_;

        //- Species fields the model loaded from disk (those not already
        //  registered with the mesh by an external solver). Held here for
        //  lifetime; access is via the mesh registry in correct().
        PtrList<volScalarField> ownedSpeciesFields_;

        //- Absorption coefficients (one per species per band)
        DynamicList<DynamicList<scalar, 1, 0> > aCoeffs_;

        //- Switch of scattering model
        Switch scattering_;

        //- Number of scattering species
        label nScattering_;

        //- Names of scattering species
        DynamicList<word, 1, 0> sSpecies_;

        //- Scattering coefficients (one per species per band)
        DynamicList<DynamicList<scalar, 1, 0> > sCoeffs_;
       
public:

    //- Runtime type information
    TypeName("linearSpecies");


    // Constructors

        //- Construct from components
        linearSpeciesExtinction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    // Destructor
    virtual ~linearSpeciesExtinction();

    // Member Operators

        //- Correct the extinction coefficients
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace optical
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
