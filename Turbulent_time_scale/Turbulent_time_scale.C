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

#include "Turbulent_time_scale.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "IncompressibleTurbulenceModel.H"
#include "RASModel.H"


namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Turbulent_time_scale, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Turbulent_time_scale,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::Turbulent_time_scale::calc()
{	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");
		
		const volScalarField& kcopy = Ucopy.db().lookupObject<volScalarField>("k");


		const volScalarField& epsiloncopy = Ucopy.db().lookupObject<volScalarField>("epsilon");

		const tmp<volScalarField> tempA(kcopy/epsiloncopy);

		const volScalarField& finalResult = tempA;

                return store(resultName_, finalResult * 1 );
}

	


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Turbulent_time_scale::Turbulent_time_scale
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

Foam::functionObjects::Turbulent_time_scale::~Turbulent_time_scale()
{}


// ************************************************************************* //
