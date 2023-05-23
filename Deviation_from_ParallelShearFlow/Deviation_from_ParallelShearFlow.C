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

#include "Deviation_from_ParallelShearFlow.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(Deviation_from_ParallelShearFlow, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        Deviation_from_ParallelShearFlow,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
//

bool Foam::functionObjects::Deviation_from_ParallelShearFlow::calc()
{	
		const volVectorField& Ucopy = lookupObject<volVectorField>("U");
		
		const tmp<volTensorField> tgradU(fvc::grad(Ucopy));	

		const volTensorField& gradU = tgradU();

		volScalarField a(IOobject("a", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedScalar("a", dimensionSet(0,2,-3,0,0,0,0), Foam::scalar(0.0)));

		forAll(a.ref(), cellI)
		{
			for(int i = 0; i < 3; i++)
			{
				for(int j = 0; j < 3; j++)
				{
					a.ref()[cellI] += (Ucopy.internalField()[cellI].component(i) * Ucopy.internalField()[cellI].component(j) * gradU.internalField()[cellI].component((3 *i) + j)); 
				}
			}
			a.ref()[cellI] = mag(a.ref()[cellI]);
		}
		

		forAll(a.boundaryFieldRef(), patchI)
                {
                        forAll(a.boundaryFieldRef()[patchI], faceI)
                        {
                                for(int i = 0; i < 3; i++)
                                {
                                        for(int j = 0; j < 3; j++)
                                        {
                                                a.boundaryFieldRef()[patchI][faceI] += (a.boundaryFieldRef()[patchI][faceI] + (Ucopy.boundaryField()[patchI][faceI].component(i) * Ucopy.boundaryField()[patchI][faceI].component(j) * gradU.boundaryField()[patchI][faceI].component((3 *i) + j)));
                                        }
                                }
                                a.boundaryFieldRef()[patchI][faceI] = mag(a.boundaryFieldRef()[patchI][faceI]);
                        }
                }
		
		return store(resultName_, a * 1); 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::Deviation_from_ParallelShearFlow::Deviation_from_ParallelShearFlow
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

Foam::functionObjects::Deviation_from_ParallelShearFlow::~Deviation_from_ParallelShearFlow()
{}


// ************************************************************************* //
