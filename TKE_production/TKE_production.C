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

#include "TKE_production.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceFields.H"
#include "fvcGrad.H"

#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "IncompressibleTurbulenceModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TKE_production, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        TKE_production,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::TKE_production::calc()
{	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");
		
		const tmp<volTensorField> tgradU(fvc::grad(Ucopy));

		const volTensorField& gradU = tgradU;

		volSymmTensorField S(symm(gradU));
		
                typedef incompressible::turbulenceModel icoModel;

		tmp<volSymmTensorField> Rcopy;

                if(mesh_.foundObject<icoModel>(turbulenceModel::propertiesName))
                {
                        const icoModel& model = mesh_.lookupObject<icoModel>(turbulenceModel::propertiesName);
                        Rcopy = model.R();
                }
                else
                {
                        FatalErrorInFunction << "Unable to find turbulence model in the database" << exit(FatalError);
                }

		volScalarField result = mag(Rcopy & S); 

		return store(resultName_, result * 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TKE_production::TKE_production
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

Foam::functionObjects::TKE_production::~TKE_production()
{}


// ************************************************************************* //
