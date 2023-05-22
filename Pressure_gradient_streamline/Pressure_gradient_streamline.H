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

Class
    Foam::functionObjects::Pressure_gradient_streamline

Description
    Calculates and outputs the second invariant of the velocity gradient tensor
    [1/s^2].

    \f[
        Pressure_gradient_streamline = 0.5(sqr(tr(\nabla U)) - tr(((\nabla U) \cdot (\nabla U))))
    \f]

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    Pressure_gradient_streamline.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Pressure_gradient_streamline_H
#define functionObjects_Pressure_gradient_streamline_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class Pressure_gradient_streamline Declaration
\*---------------------------------------------------------------------------*/

class Pressure_gradient_streamline
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the Pressure_gradient_streamline field and return true if successful
        virtual bool calc();


//	volVectorField uMean(IOobject("uMean", Ucopy.mesh().time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh_);// dimensionedVector("uMean", dimensionSet(0,1,-1,0,0,0,0), Foam::vector(0,0,0)));

public:

    //- Runtime type information
    TypeName("Pressure_gradient_streamline");
    
 
    static const word modelName;

//     const volVectorField& Ucopy = lookupObject<volVectorField>("U");

//     volVectorField uMean = Ucopy;


    // Constructors

        //- Construct from Time and dictionary
        Pressure_gradient_streamline
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Pressure_gradient_streamline();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
