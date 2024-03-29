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
    Foam::functionObjects::contactAngleEvaluation

SourceFiles
    contactAngleEvaluation.C

Description
    Calculate the contact angle in 3D simulation domain
    
    From the interface cells at the wall,
    identifies the contact line cell by reconstructing the signed distances at the vertices of the cell
    and calculates the distance of the contact line from the central axis of the drop.

    
Author
    Muhammad Hassan Asghar
    asghar@mma.tu-darmstadt.de, hassan.asghar@tu-darmstadt.de
    Mathematical Modeling and Analysis Group
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#ifndef contactAngleEvaluation_H
#define contactAngleEvaluation_H

#include "functionObject.H"

#include "volFields.H"

#include "OFstream.H"


 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    /*---------------------------------------------------------------------------*\
                       Class contactAngleEvluation Declaration
    \*---------------------------------------------------------------------------*/

    class contactAngleEvaluation
        :
        public functionObject {
            // Private Data

            //-Time
            const Time & time_;

            //- Mesh const reference.
            const fvMesh & mesh_;

            //- Name of the phase indicator field.
            word fieldName_;
            
            point origin_;

            //- Reference to the phase indicator.
            const volScalarField & alpha1_;

            //Reference to PLIC normal
            const volVectorField & plicNormals_;
            //Reference to PLIC centre
            const volVectorField & plicCentres_;

            // stores the contact angle field
            volScalarField contactAngles_;

            // stores the contact radius field
            volScalarField rf_;

            // stores the wisps as field for visualization
            volScalarField wisp_;
            
            //min contact angle in the domain
            scalar contactAnglemin_;

            // maximum contact angle in the domain
            scalar contactAnglemax_;

            //Maximum contact radius
            scalar rfMax_;
          
	               
            //DynamicList <scalar> contactAngle;
            
            //tolerance for cell to be the interface cell
            scalar alphaTol_ = 1e-8;
            
            // wisp tolerance to mark if a cellis a wisp or not
            scalar wispTol_ = 1e-03;

	        //- Visualization 
            //mutable volScalarField contactLineAngle_;

            public:

                //- Runtime type information
                TypeName("contactAngleEvaluation");

            // Constructors

            //- Construct from Time and dictionary
            contactAngleEvaluation
                (
                    const word & name,
                        const Time & time,
                            const dictionary & dict
                );

            //- No copy construct
            contactAngleEvaluation(const contactAngleEvaluation & ) = delete;

            //- No copy assignment
            void operator = (const contactAngleEvaluation & ) = delete;

            //- Destructor
            virtual~contactAngleEvaluation() =
                default;

            // Member Functions
            //- Access time
            const Time & getTime() const {
                return time_;
            }

            //- Read the contactAngleEvaluation data
            virtual bool read(const dictionary & dict);

            //- Execute, currently does nothing
            virtual bool execute();

            //- Called at the start of the time-loop
            virtual bool start();

       	    void maxRf (const volScalarField& rf_);

	     //- Calculate the contact angle
            void calContactAngle();

	    // marks the bulk cells that are wisps (Maric et. al.)
            void checkWisps();

            //- Calculate the rf_ //Distance from central axis of the drop to the contact line cell
            void calRf(label bCell, const volVectorField& C, int &rcount);

            //- Determine whether a face contains a contact line
            bool hasContactLine
            (
                label faceI,
                const fvPatch& patch
            ) const;

           // void isAWisp (scalar threshold);

            void report();
            //- Execute at the final time-loop, currently does nothing
            virtual bool end();

            //- Write the contactAngleEvaluation
            virtual bool write();
            
        };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
