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

#include "deltaEddingtonModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace photoBio
    {
        defineTypeNameAndDebug(deltaEddingtonModel, 0);

        addToRunTimeSelectionTable
        (
            inScatterModel,
            deltaEddingtonModel,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::photoBio::deltaEddingtonModel::deltaEddingtonModel
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

     coeffsDict_.lookup("nBand") >> nBands_;

     phaseFuncCoef_.setSize(nBands_);

     coeffsDict_.lookup("phaseFuncCoef") >> phaseFuncCoef_;

     coeffsDict_.lookup("forwardScatterFactor") >> forwardScatterFactor_;
      

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::photoBio::deltaEddingtonModel::~deltaEddingtonModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


 Foam::scalar  Foam::photoBio::deltaEddingtonModel::correct
(
          const  scalar  rayCos,
          const   label  iBand
) const
{
	   scalar y ; 
       if( rayCos == 1 )   y = 2*forwardScatterFactor_[iBand]+ (1-forwardScatterFactor_[iBand])*(1 +  phaseFuncCoef_[iBand]);
       else            y = (1 - forwardScatterFactor_[iBand]) * (1 +   phaseFuncCoef_[iBand]*rayCos );
       
	   return y;

}





// ************************************************************************* //
