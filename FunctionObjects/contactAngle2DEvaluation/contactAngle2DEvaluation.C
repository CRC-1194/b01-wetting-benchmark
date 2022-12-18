/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "contactAngle2DEvaluation.H"

#include "addToRunTimeSelectionTable.H"

#include <iostream>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {

    defineTypeNameAndDebug(contactAngle2DEvaluation, 0);
    addToRunTimeSelectionTable(functionObject, contactAngle2DEvaluation, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    contactAngle2DEvaluation::contactAngle2DEvaluation(
            const word & name,
            const Time & time,
            const dictionary & dict
        ):
        functionObject(name),
        time_(time),
        mesh_(time.lookupObject < fvMesh > (polyMesh::defaultRegion)),
        fieldName_(dict.lookup("phaseIndicator")),
        patchName_(dict.lookup("patchName")),
        alpha1_
        (
            mesh_.lookupObject < volScalarField >
            (
                fieldName_
            )
        ),
        plicNormals_
        (
            mesh_.lookupObject < volVectorField >
            (
               IOobject::groupName("interfaceNormal", alpha1_.group())
            )
        ),
        plicCentres_
        (
            mesh_.lookupObject < volVectorField >
            (
                IOobject::groupName("interfaceCentre", alpha1_.group())
            )
        ),
        contactAngles_
        (
            IOobject
            (
                "contactAngles_",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "contactAngles_",
                dimless,
                0.0
            )
        ),
        outputFile("contactAngles.csv"),
        contactAnglemin_(0), 
        contactAnglemax_(0)
    {
        outputFile << "TIME,CONTACT_ANGLE_LEFT,CONTACT_ANGLE_RIGHT,CONTACT_ANGLE_ADVANCING_DEG,CONTACT_ANGLE_RECEDING_DEG\n";
        read(dict);
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool contactAngle2DEvaluation::read(const dictionary & dict) {
        //notImplemented("contactAngle2DEvaluation::read(const dictionary&)");
        //return execute();
        calContactAngle();
        contactAngles_.write();
        return true;
    }

    bool contactAngle2DEvaluation::execute() {
        // Calculate advancing and receding contact angles at time t > 0.
        calContactAngle();
//        report();
        return true;
    }

    //From all interface cell
    // Identifies the contact line cell at the wall
    // Calculate the maximum and minimum contact angle
    void contactAngle2DEvaluation::calContactAngle() {
        contactAngles_*= scalar(0);
        //const auto& faceOwner = mesh_.faceOwner();
        double patchCellCentreXComponent, patchCellCentreXComponent_temp = 0.0;

        //get the patch ID number
        label patchID = mesh_.boundaryMesh().findPatchID(patchName_);

        const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
        const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();
        
        //- Get boundary mesh
        const fvBoundaryMesh & boundary = mesh_.boundary();
        //- Get desired path mesh from boundary mesh
        const fvPatch & patch = boundary[patchID];
        // Get normals to the patch cell face
        const vectorField & nHatf = patch.nf();

        const vectorField & nC = patch.Cf();
        int count =0;
        int reqCellID_min = 0;
        int reqCellID_max = 0;

        // Detects the cells with contact line
        // Strictly applied to 2D only
        forAll (plicNormalPatchInternalField, faceI)
        {
          if (plicNormalPatchInternalField[faceI] != vector::zero){
            if(count==0)
            {
                patchCellCentreXComponent_temp = nC[faceI].component(vector::X);
                patchCellCentreXComponent = patchCellCentreXComponent_temp;
                count++;
                reqCellID_max = faceI;
                reqCellID_min = faceI;
            }
            else
            {
                patchCellCentreXComponent_temp = nC[faceI].component(vector::X);
                if(patchCellCentreXComponent < patchCellCentreXComponent_temp)
                {
                    patchCellCentreXComponent = patchCellCentreXComponent_temp;
                    reqCellID_min = faceI;
                }
                else if(patchCellCentreXComponent > patchCellCentreXComponent_temp)
                {
                    patchCellCentreXComponent = patchCellCentreXComponent_temp;
                    reqCellID_max = faceI;
                }
            }  
          }  
           
        }

        // contact angle (maximum and minimum)from the interface unit normal and boundary outer unit normal vector
        contactAnglemin_ = acos((plicNormalPatchInternalField[reqCellID_min]) & nHatf[reqCellID_min] / (mag (plicNormalPatchInternalField[reqCellID_min]) * mag(nHatf[reqCellID_min])) ) * 180.0 / M_PI; 
        contactAngles_[mesh_.boundary()[patchID].faceCells()[reqCellID_min]] = contactAnglemin_;
        contactAnglemax_ = acos((plicNormalPatchInternalField[reqCellID_max]) & nHatf[reqCellID_max] / (mag (plicNormalPatchInternalField[reqCellID_max]) * mag(nHatf[reqCellID_max])) ) * 180.0 / M_PI; 
        contactAngles_[mesh_.boundary()[patchID].faceCells()[reqCellID_max]] = contactAnglemax_;
        Info << contactAnglemin_ << " "  << contactAnglemax_;
        outputFile <<  time_.timeOutputValue() << "," // physical time 
                   << contactAnglemax_ << "," << contactAnglemin_ << endl;
        
    }

    //- Report at the console output / log file
    void contactAngle2DEvaluation::report() { 
        
        //Info << endl << endl << "Contact angle is " << contactAngle_ <<
          //  " at time step " << time_.value() << endl << endl;
    }

    bool contactAngle2DEvaluation::start() {
        //execute();
        return true;
    }

    bool contactAngle2DEvaluation::end() {
        return true;
    }

    bool contactAngle2DEvaluation::write() {  
        if (time_.writeTime()) {
            // Create the output path directory
            contactAngles_.write();
            reduce(contactAnglemin_, minOp<scalar>());
            reduce(contactAnglemax_, maxOp<scalar>());
            fileName outputDir = "postProcessing";
            // Createe the directory
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr < OFstream > outputFilePtrMin;
            autoPtr < OFstream > outputFilePtrMax;
            // Open the file in the newly created directory
            outputFilePtrMin.reset(new OFstream(outputDir/"contactangleMin"+".csv",
                                IOstreamOption(),
                                true));
            outputFilePtrMax.reset(new OFstream(outputDir/"contactangleMax"+".csv",
                                IOstreamOption(),
                                true));

            // Write stuff
            outputFilePtrMin() << time_.value()<< "," << contactAnglemin_<< endl;
            outputFilePtrMax() << time_.value()<< "," << contactAnglemax_<< endl;

        }   
      return true;
    }

}
// ************************************************************************* //
