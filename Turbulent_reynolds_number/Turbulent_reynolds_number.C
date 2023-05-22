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

#include "Turbulent_reynolds_number.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Turbulent_reynolds_number, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Turbulent_reynolds_number,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::Turbulent_reynolds_number::calc()
{	
	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");	

		const volScalarField& kcopy = Ucopy.db().lookupObject<volScalarField>("k");

			
//	        const volScalarField& dist =  Foam::wallDist(mesh_).y(); 		

		const volScalarField& nucopy = Ucopy.db().lookupObject<volScalarField>("nu");


		return store(resultName_, (sqrt(kcopy) * Foam::wallDist(mesh_).y()) / (50 * nucopy));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Turbulent_reynolds_number::Turbulent_reynolds_number
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

Foam::functionObjects::Turbulent_reynolds_number::~Turbulent_reynolds_number()
{}


// ************************************************************************* //
