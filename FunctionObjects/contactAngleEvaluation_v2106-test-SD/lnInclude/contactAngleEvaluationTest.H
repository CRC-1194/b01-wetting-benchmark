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

Class
    Foam::functionObjects::wettedArea

Group
    grpFieldFunctionObjects

Description
    <minimal description of the function object>

    <equation>
    \f[
        x = x_{ref}^x + \rho \omega
    \f]

    <variable-explanation table>
    where
    \vartable
      \rho        | <explanation> [units, e.g. kg/m3]
      \omega      | \f$ \nabla \cdot \vec U \f$
      ...         | ...
    \endvartable

    <inline equation>
    where \f$ x_k \f$ is ...

    <input-output table>
    \table
      Operand      | Type                 | Location
      input        | {vol,surface}\<Type\>Field(s) <!--
               --> |$FOAM_CASE/\<time\>/\<inpField\>s
      output file  | dat                  <!--
              -->  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field | volScalarField       | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wettedArea1
    {
        // Mandatory entries (unmodifiable)
        type           wettedArea;
        libs           (wettedAreaFunctionObject);

        // Mandatory entries (runtime modifiable)
        ...

        // Mandatory (inherited) entries (unmodifiable)
        ...

        // Mandatory (inherited) entries (runtime unmodifiable)
        ...

        // Optional entries (unmodifiable)
        ...

        // Optional entries (runtime modifiable)
        boolData       <bool>;
        labelData      <label>;
        wordData       <word>;
        scalarData     <scalar>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                      | Type | Req'd | Dflt
      type       | Type name: wettedArea        | word |  yes  | -
      libs       | Library name: wettedAreaFunctionObject  <!--
             --> | word  | yes | -
      boolData   | <explanation>                    | bool |  yes  | -
      labelData  | <explanation>                    | label |  yes | -
      wordData   | <explanation>                    | word |  yes  | -
      scalarData | <explanation>                    | scalar |  no | 1.0
      wordListData  | <explanation>                 | wordList | yes | -
    \endtable

    Options for the \c ENTRY entry:
    \verbatim
        <option1>
        <option2> | <explanation>
      ...
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link fieldsExpression.H \endlink
     - \link writeFile.H \endlink
     ...

    <if \c postProcess is applicable>
    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func wettedArea
    \endverbatim

    <if \c postProcess is not applicable>
    Usage by the \c postProcess utility is not available.

Note
    - <note1>
    - <note2>
    ...

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::wettedArea
    ...

SourceFiles
    contactAngleEvaluationTest.C
    wettedAreaTEMPLATES.C
    ...

\*---------------------------------------------------------------------------*/

#ifndef contactAngleEvaluationTest_H
#define contactAngleEvaluationTest_H

#include "functionObject.H"

#include "volFields.H"

#include "OFstream.H"


 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

    /*---------------------------------------------------------------------------*\
                       Class contactAngleEvluation Declaration
    \*---------------------------------------------------------------------------*/

    class contactAngleEvaluationTest
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

            const volVectorField & plicNormals_;

            const volVectorField & plicCentres_;

            //- Wetted Area total.
            
            volScalarField contactAngles_;

            volScalarField rf_;

            volScalarField wisp_;
            
            scalar contactAnglemin_;

            scalar contactAnglemax_;

            scalar rfMax_;
          

            DynamicList <scalar> contactAngle;
            //tolerance for cell to be the interface cell
            scalar alphaTol_ = 1e-8;

            scalar wispTol_ = 1e-03;

            public:

                //- Runtime type information
                TypeName("contactAngleEvaluationTest");

            // Constructors

            //- Construct from Time and dictionary
            contactAngleEvaluationTest
                (
                    const word & name,
                        const Time & time,
                            const dictionary & dict
                );

            //- No copy construct
            contactAngleEvaluationTest(const contactAngleEvaluationTest & ) = delete;

            //- No copy assignment
            void operator = (const contactAngleEvaluationTest & ) = delete;

            //- Destructor
            virtual~contactAngleEvaluationTest() =
                default;

            // Member Functions
            //- Access time
            const Time & getTime() const {
                return time_;
            }

            //- Read the contactAngleEvaluationTest data
            virtual bool read(const dictionary & dict);

            //- Execute, currently does nothing
            virtual bool execute();

            //- Called at the start of the time-loop
            virtual bool start();

       	    void maxRf (const volScalarField& rf_);

	     //- Calculate the contact angle
            void calContactAngle();

            void checkWisps();

            //- Calculate the rf_ //Distance from central axis of the drop to the contact line cell
            void calRf(label bCell, const volVectorField& C, int &rcount);

            void report();
            //- Execute at the final time-loop, currently does nothing
            virtual bool end();

            //- Write the contactAngleEvaluationTest
            virtual bool write();
            
        };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
