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

\*---------------------------------------------------------------------------*/

#include "anisotropicScattering.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace optical
    {
        defineTypeNameAndDebug(anisotropicScattering, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
            anisotropicScattering,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optical::anisotropicScattering::anisotropicScattering
(
    const dictionary& dict
  //  , const fvMesh& mesh
)
:
 //   inScatterModel(dict, mesh),
     inScatterModel(dict),
      coeffsDict_(dict.subDict(typeName + "Coeffs"))
 //  , totalWaveLength_(0)
{
	
//    label iBand = 0;
 //   const dictionary& functionDicts = dict.subDict(typeName +"Coeffs");
    
     coeffsDict_.lookup("nBand") >> nBands_;

     phaseFuncCoef_.setSize(nBands_);

     coeffsDict_.lookup("phaseFuncCoef") >> phaseFuncCoef_;
    
    /*          
    forAllConstIter(dictionary, functionDicts, iter)
    {
        // safety:
        if (!iter().isDict())
        {
            continue;
        }
        phaseFunc_[iBand].setCapacity(phaseFunctionAngleNum_);   
        const dictionary& dict = iter().dict();
        dict.lookup("phaseFunction") >> phaseFunc_[iBand];
        iBand++;
    }
*/

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optical::anisotropicScattering::~anisotropicScattering()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


 Foam::scalar  Foam::optical::anisotropicScattering::correct
(
          const  scalar  rayCos,
          const   label  iBand
) const
{

       scalar y = 1 +  phaseFuncCoef_[iBand]*rayCos;
       
	   return y;

}





// ************************************************************************* //
