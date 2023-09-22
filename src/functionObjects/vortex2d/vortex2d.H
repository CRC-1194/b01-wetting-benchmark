/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::vortex2d

Description
    Provides options to set the velocity and flux fields as a function of time.

    Useful for testing advection behaviour of numerical schemes by e.g.
    imposing solid body rotation, vortex flows.  

SourceFiles
    vortex2d.C
    
Author
    Muhammad Hassan Asghar
    asghar@mma.tu-darmstadt.de, hassan.asghar@tu-darmstadt.de
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_vortex2d_H
#define functionObjects_vortex2d_H

#include "fvMeshFunctionObject.H"
#include "Function1.H"
#include "Enum.H"
#include "point.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class vortex2d Declaration
\*---------------------------------------------------------------------------*/

class vortex2d
:
    public fvMeshFunctionObject
{
    // Private Enumerations

        //- Options for the operating mode
        enum class modeType
        {
            MYVORTEX2D
        };

        //- Names for modeType
        static const Enum<modeType> modeTypeNames;


   // Private Data

        //- Operating mode
        modeType mode_;

        //- Name of velocity field
        word UName_;

        //- Name of density field
        word rhoName_;

        //- Name of flux field
        word phiName_;

        //- Reverse time
        scalar reverseTime_;

        //- Scaling function
        autoPtr<Function1<scalar>> scalePtr_;


        // Rotation

            //- Origin
            point origin_;

            //- Rotation tensor for rotational mode
            tensor R_;

            //- Rotational speed function
            autoPtr<Function1<scalar>> omegaPtr_;


        // Function

            //- Velocity function
            autoPtr<Function1<vector>> velocityPtr_;


    // Private Member Functions

        //- Set the flux field
        void setPhi(const volVectorField& U);


public:

    //- Runtime type information
    TypeName("vortex2d");


    // Constructors

        //- Construct from Time and dictionary
        vortex2d
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        vortex2d(const vortex2d&) = delete;

        //- No copy assignment
        void operator=(const vortex2d&) = delete;


    //- Destructor
    virtual ~vortex2d() = default;


    // Member Functions

        //- Read the vortex2d data
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Calculate the vortex2d and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
