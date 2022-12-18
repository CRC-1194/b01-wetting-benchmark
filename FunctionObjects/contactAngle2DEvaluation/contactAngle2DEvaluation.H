/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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
    Foam::functionObjects::contactAngle2DEvaluation

SourceFiles
    contactAngle2DEvaluation.C
    
Description
    Calculate the contact angle in 2D simulations
    
    From the interface cells at the wall,
    identifies the contact line cell and writes
    the maximum and minimum contact angle
    
Author
    Muhammad Hassan Asghar
    asghar@mma.tu-darmstadt.de, hassan.asghar@tu-darmstadt.de
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

    ...

\*---------------------------------------------------------------------------*/

#ifndef contactAngle2DEvaluation_H
#define contactAngle2DEvaluation_H

#include "functionObject.H"

#include "volFields.H"

#include "OFstream.H"


 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    /*---------------------------------------------------------------------------*\
                       Class contactAngleEvluation Declaration
    \*---------------------------------------------------------------------------*/

    class contactAngle2DEvaluation
        :
        public functionObject {
            // Private Data

            //-Time
            const Time & time_;

            //- Mesh const reference.
            const fvMesh & mesh_;

            //- Name of the phase indicator field.
            word fieldName_;
            
            //- Name of the patch
            word patchName_;

            //- Reference to the phase indicator.
            const volScalarField & alpha1_;

	    // Normals to the PLIC plane
            const volVectorField & plicNormals_;

            // Centres of the PLIC component
            const volVectorField & plicCentres_;

            //- Contact Angle scalar field for visualization
            volScalarField contactAngles_;

	     // Output file stream for the advancing and receding contact angles. 
            OFstream outputFile;  
            
            //Minimum contact angle in the domain
            scalar contactAnglemin_;

	    //Maximum contact angle in the domain
            scalar contactAnglemax_;

            DynamicList <scalar> contactAngle;
            
            //tolerance for cell to be the interface cell
            scalar alphaTol_ = 1e-8;

            public:

                //- Runtime type information
                TypeName("contactAngle2DEvaluation");

            // Constructors

            //- Construct from Time and dictionary
            contactAngle2DEvaluation
                (
                    const word & name,
                        const Time & time,
                            const dictionary & dict
                );

            //- No copy construct
            contactAngle2DEvaluation(const contactAngle2DEvaluation & ) = delete;

            //- No copy assignment
            void operator = (const contactAngle2DEvaluation & ) = delete;

            //- Destructor
            virtual~contactAngle2DEvaluation() =
                default;

            // Member Functions
            //- Access time
            const Time & getTime() const {
                return time_;
            }

            //- Read the contactAngle2DEvaluation data
            virtual bool read(const dictionary & dict);

            //- Execute, currently does nothing
            virtual bool execute();

            //- Called at the start of the time-loop
            virtual bool start();

            //- Calculate the wetted area
            void calContactAngle();

            void report();
            //- Execute at the final time-loop, currently does nothing
            virtual bool end();

            //- Write the contactAngle2DEvaluation
            virtual bool write();
            
        };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
