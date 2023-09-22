/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR, AFFILIATION
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

SourceFiles
    updateFlux.C
    
Description:
	Calculates the wetted area of the specified patch 
	and write it in a .csv file

Author
    Muhammad Hassan Asghar
    asghar@mma.tu-darmstadt.de, hassan.asghar@tu-darmstadt.de
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#ifndef updateFlux_H
#define updateFlux_H

#include "functionObject.H"

#include "volFields.H"

#include "OFstream.H"
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    /*---------------------------------------------------------------------------*\
                       Class updateFlux Declaration
    \*---------------------------------------------------------------------------*/

    class updateFlux
        :
        public functionObject {
            // Private Data

            //-Time
            const Time & time_;

            //- Mesh const reference.
            const fvMesh & mesh_;

            //- Name of the phase indicator field.
            word fieldName_;

            //Acceleration of the plate retraction
            vector acceleration_;

            //Initial height of the cells :  below this height, update of U and flux field will be done at each time step
            // This will act as an imaginary wall which is retracted downwards with a velocity given by the acceleration_
            // 
            vector initialHeight_;

            //- Reference to the phase indicator.
            const volScalarField & alpha1_;

            //- Reference to the cells whose flux is being updated .
            volScalarField updatedCells_;

            //- updated flux vlaue total.
            scalar updateFlux_;

            public:

                //- Runtime type information
                TypeName("updateFlux");

            // Constructors

            //- Construct from Time and dictionary
            updateFlux
                (
                    const word & name,
                        const Time & time,
                            const dictionary & dict
                );

            //- No copy construct
            updateFlux(const updateFlux & ) = delete;

            //- No copy assignment
            void operator = (const updateFlux & ) = delete;

            //- Destructor
            virtual~updateFlux() =
                default;

            // Member Functions
            //- Access time
            const Time & getTime() const {
                return time_;
            }

            //- Read the updateFlux data
            virtual bool read(const dictionary & dict);

            //- Execute, currently does nothing
            virtual bool execute();

            //- Called at the start of the time-loop
            virtual bool start();

            //- Calculate the wetted area
            void calculateupdateFlux();

            void report();
            //- Execute at the final time-loop, currently does nothing
            virtual bool end();

            //- Write the updateFlux
            virtual bool write();
        };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
