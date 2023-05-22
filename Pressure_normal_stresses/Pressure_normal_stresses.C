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

#include "Pressure_normal_stresses.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Pressure_normal_stresses, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Pressure_normal_stresses,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::Pressure_normal_stresses::calc()
{	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");
		
		const volScalarField& pcopy = lookupObject<volScalarField>("p");

		const tmp<volVectorField> tgradP(fvc::grad(pcopy));	

		const volVectorField& gradP = tgradP();

		const tmp<volScalarField> temp(gradP & gradP);

		const volScalarField& tempa = temp;

		const dimensionedScalar& rhocopy = Ucopy.db().lookupObject<IOdictionary>("transportProperties").lookup("rho");

		const tmp<volScalarField> finalResult(sqrt(tempa)); // (rho * rho) removed, why was it there?

		return store(resultName_, finalResult * 1); 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Pressure_normal_stresses::Pressure_normal_stresses
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

Foam::functionObjects::Pressure_normal_stresses::~Pressure_normal_stresses()
{}


// ************************************************************************* //
