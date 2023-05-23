/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "Streamline_curvature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "ddt.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Streamline_curvature, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Streamline_curvature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::Streamline_curvature::calc()
{	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");
		
		const surfaceScalarField& phicopy = Ucopy.db().lookupObject<surfaceScalarField>("phi");

		volScalarField a = mag(fvc::ddt(Ucopy) + (fvc::div(phicopy, Ucopy))); 	

		forAll(Ucopy.internalField(), faceI)
		{
		         a.ref()[faceI] = a.ref()[faceI] / magSqr(Ucopy.internalField()[faceI]); 
		}

		forAll(Ucopy.boundaryField(), patchI)
		{
		         forAll(Ucopy.boundaryField()[patchI], faceI)
		         {
			       if ((Ucopy.boundaryField()[patchI][faceI].x()) != 0 || (Ucopy.boundaryField()[patchI][faceI].y()) != 0 || (Ucopy.boundaryField()[patchI][faceI].z()) != 0)
		               { 																									                a.boundaryFieldRef()[patchI][faceI] = (a.boundaryFieldRef()[patchI][faceI] / (magSqr(Ucopy.boundaryField()[patchI][faceI])));	
	  		       }
	                       else
	                       {
	                           a.boundaryFieldRef()[patchI][faceI] = 0.0;
			       }
			 }
		}


		return store(resultName_, a * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Streamline_curvature::Streamline_curvature
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName(typeName, "U");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::Streamline_curvature::~Streamline_curvature()
{}


// ************************************************************************* //
