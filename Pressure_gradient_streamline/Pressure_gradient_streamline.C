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

#include "Pressure_gradient_streamline.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Pressure_gradient_streamline, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Pressure_gradient_streamline,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::Pressure_gradient_streamline::calc()
{	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");
		
		const volScalarField& pcopy = lookupObject<volScalarField>("p");

		const tmp<volVectorField> tgradP(fvc::grad(pcopy));	

		const volVectorField& gradP = tgradP();

		const tmp<volScalarField> temp(Ucopy & gradP);

		const volScalarField& finalResult = temp;
               
		return store(resultName_, finalResult * 1 ); 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Pressure_gradient_streamline::Pressure_gradient_streamline
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

Foam::functionObjects::Pressure_gradient_streamline::~Pressure_gradient_streamline()
{}


// ************************************************************************* //
