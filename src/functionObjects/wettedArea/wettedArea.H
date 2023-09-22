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
    wettedArea.C
    
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

#ifndef wettedArea_H
#define wettedArea_H

#include "functionObject.H"

#include "volFields.H"

#include "OFstream.H"
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    /*---------------------------------------------------------------------------*\
                       Class wettedArea Declaration
    \*---------------------------------------------------------------------------*/

    class wettedArea
        :
        public functionObject {
            // Private Data

            //-Time
            const Time & time_;

            //- Mesh const reference.
            const fvMesh & mesh_;

            //- Name of the phase indicator field.
            word fieldName_;

            //- Reference to the phase indicator.
            const volScalarField & alpha1_;

            //- Wetted Area total.
            scalar wettedArea_;

            public:

                //- Runtime type information
                TypeName("wettedArea");

            // Constructors

            //- Construct from Time and dictionary
            wettedArea
                (
                    const word & name,
                        const Time & time,
                            const dictionary & dict
                );

            //- No copy construct
            wettedArea(const wettedArea & ) = delete;

            //- No copy assignment
            void operator = (const wettedArea & ) = delete;

            //- Destructor
            virtual~wettedArea() =
                default;

            // Member Functions
            //- Access time
            const Time & getTime() const {
                return time_;
            }

            //- Read the wettedArea data
            virtual bool read(const dictionary & dict);

            //- Execute, currently does nothing
            virtual bool execute();

            //- Called at the start of the time-loop
            virtual bool start();

            //- Calculate the wetted area
            void calculateWettedArea();

            void report();
            //- Execute at the final time-loop, currently does nothing
            virtual bool end();

            //- Write the wettedArea
            virtual bool write();
        };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
