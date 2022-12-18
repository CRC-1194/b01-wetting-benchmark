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

#include "contactAngleEvaluationTest.H"

#include "addToRunTimeSelectionTable.H"

#include <iostream>

#include "vectorList.H"

#include <cmath>




// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {

    defineTypeNameAndDebug(contactAngleEvaluationTest, 0);
    addToRunTimeSelectionTable(functionObject, contactAngleEvaluationTest, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    contactAngleEvaluationTest::contactAngleEvaluationTest(
            const word & name,
            const Time & time,
            const dictionary & dict
        ):
        functionObject(name),
        time_(time),
        mesh_(time.lookupObject < fvMesh > (polyMesh::defaultRegion)),
        fieldName_(dict.lookup("phaseIndicator")),
        origin_(dict.get<point>("centre")),
        alpha1_(
            mesh_.lookupObject < volScalarField >
            (
                fieldName_
            )
        ),
        plicNormals_(
            mesh_.lookupObject < volVectorField >
            (
               IOobject::groupName("interfaceNormal", alpha1_.group())
            )
        ),
        plicCentres_(
            mesh_.lookupObject < volVectorField >
            (
                IOobject::groupName("interfaceCentre", alpha1_.group())
            )
        ),
        contactAngles_
        (
            IOobject
            (
                "contactAngles",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimless,
                0.0
            )
),
        rf_
        (
            IOobject
            (
                "rf",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimLength,
                0.0
            )
        ),
        wisp_
        (
            IOobject
            (
                "wisp",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimless,
                0.0
            )
        ),
        contactAnglemin_(0),
        contactAnglemax_(0),
        rfMax_(0)
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


    bool contactAngleEvaluationTest::read(const dictionary & dict) {
        notImplemented("contactAngleEvaluationTest::read(const dictionary&)");
        return execute();
    }

    bool contactAngleEvaluationTest::execute() {
        checkWisps();
        calContactAngle();
        report();
        return true;
    }

    // It calculates the maximum contact radius from all the contact points in 3D domain
    void contactAngleEvaluationTest::maxRf (const volScalarField& rf_)
    {
        rfMax_ = rf_[0];
        forAll(rf_, r)
        {
                if(rf_[r] > rfMax_)
                {
                        rfMax_ = rf_[r];
                }
        }
    }

    // Calculates the contact point distance from the central axis -> called contact radius and store it as a vol scalar field
    void contactAngleEvaluationTest::calRf(label bCell, const volVectorField& C, int &rcount)
    {
        // provides the distance between the central axis and the contact line cell centre
        label cellID = bCell;
        point origin = origin_;
        point contactLinePoint = C[cellID ];
        point positionOnCentralAxis = {origin[0], origin[1], contactLinePoint[2]};
       
        if (mag((contactLinePoint - positionOnCentralAxis)) != 0){
            rf_[cellID] = (contactLinePoint - positionOnCentralAxis) & ((contactLinePoint - positionOnCentralAxis) / mag((contactLinePoint - positionOnCentralAxis)));
        }
    }


    //Calculates the contact angle in the interface cells, discarding the wisps
    void contactAngleEvaluationTest::calContactAngle() {
        //Mesh connectivity
        const pointField& points = mesh_.points();   // Node coordinates
        const faceList& faces = mesh_.faces();                   // Face to node
        const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion
        const volVectorField& C = mesh_.C();         // Cell center coordinates
        const auto& faceOwner = mesh_.faceOwner();
        contactAngles_*= scalar(0);
        rf_*= scalar(0);


        forAll(mesh_.boundary(), patch)
        {
            const word& patchName = mesh_.boundary()[patch].name();            // Boundary patch name

            if (patchName == "bottom" || patchName == "sphere")
            {
                label patchID = mesh_.boundaryMesh().findPatchID(patchName);
                const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
                const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();
                int count=0;
                int rcount = 0;
                rfMax_ = 0;
                 //loop over all faces of the boundary patch 
                forAll(mesh_.boundary()[patch], facei)
                {
                    // Get normals to the patch cell face
                    const vectorField& nHatf = mesh_.boundary()[patch].nf();
                    const label& face = boundaryMesh[patch].start() + facei;        // Face index
                    int bCell = faceOwner[facei + mesh_.boundary()[patch].start()]; //cell ID
                
                    // check if the cell is an interface cell
                    if(alpha1_[bCell] < (1.0-wispTol_) && alpha1_[bCell] > wispTol_ )
                    {
                        std::vector <double> signedDistancesOfPatchFace; // to store sgined distances of the interface path face
                        bool isContactLine = false;
                        int countSign = 0;
                        double contactAngle = 0.0;

                        // Signed distance  to detect the presence of the contact line. If the sign changes in even times ->  contact line
                        forAll(faces[face], vertexID) //patch face
                        {
                            point vertex = points[faces[face][vertexID]]; // vertices of a patch face with interface
                            double sd_ = plicNormalPatchInternalField[facei] & (plicCentrePatchInternalField[facei] - vertex);
                            signedDistancesOfPatchFace.push_back(sd_);
                        }
                        forAll(signedDistancesOfPatchFace, sd)
                        {
                            auto ctr = (sd + 1) % signedDistancesOfPatchFace.size();
                            if (std::signbit(signedDistancesOfPatchFace[sd]) != std::signbit(signedDistancesOfPatchFace[ctr]))
                            {
                                countSign++;
                            }
                        }
                        if(countSign!=0) // also counts for vertex intersection (2 for edge-edge intersection 4 for vertex-vertex, 2 for egde intersection), //3 for vertex-edge
                        {
                            isContactLine = true;
                            calRf(bCell, C, rcount);
                        }
                        if(isContactLine)// detects the contact line
                        { 
                            if(!wisp_[bCell]) // also checks if not a wisp 
                            {
                                if((mag (plicNormalPatchInternalField[facei]) * mag(nHatf[facei])) !=0){
                                    contactAngle = acos((plicNormalPatchInternalField[facei]) & nHatf[facei] / (mag (plicNormalPatchInternalField[facei]) * mag(nHatf[facei])) ) * 180.0 / M_PI;
                                    contactAngles_[bCell] = contactAngle;
					            }
                                        
                                if(count==0)
                                {
                                    //contactAnglemax_ = contactAngle;
                                    contactAnglemin_ = contactAngle;
                                    count++;
                                }
                                if(contactAngle < contactAnglemin_)
                                {
                                    contactAnglemin_ = contactAngle;
                                }

                            }
                        }
                    }
                }

            }
        }
    }

//Check if an interface cell is a wisp. Marks all the wisp as 1 for visualization
    void contactAngleEvaluationTest::checkWisps(){
       wisp_*= scalar(0);
        const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion 
        const auto& patches = mesh_.boundary();
        const auto& faceOwner = mesh_.faceOwner();

        forAll(patches, patch)
        {
            const word& patchName = patches[patch].name();            // Boundary patch name
            if(patchName == "bottom" || patchName == "sphere")
            {
                label patchID = boundaryMesh.findPatchID(patchName);
                const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
                const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();

                forAll(patches[patch], facei)
                {
                    const label& ownerCellId = faceOwner[facei + patches[patch].start()]; //cell ID
                    //if (plicNormalPatchInternalField[facei] != vector::zero)
                    if(alpha1_[ownerCellId] < (1.0-alphaTol_) && alpha1_[ownerCellId] > alphaTol_ )
                    {
                        const labelList& faces = boundaryMesh[patch].faceFaces()[facei];
                        bool isAWisp = true;
                        forAll(faces, face)
                        {
                            const auto& faceNId = faceOwner[faces[face] + patches[patch].start()]; //cell ID
                            if(alpha1_[faceNId] < (1.0- wispTol_) ){ //full cell
                                isAWisp = false;
                                break;
                            }
                        }
                        if(isAWisp){
                            wisp_[ownerCellId] =1.0;
                        }
                    }
                }
            }
        }
    }

    //- Report at the console output / log file
    void contactAngleEvaluationTest::report() {
        //Info << endl << endl << "Contact angle is " << contactAngle_ <<
          //  " at time step " << time_.value() << endl << endl;
    }

    bool contactAngleEvaluationTest::start() {
        //execute();
        return true;
    }

    bool contactAngleEvaluationTest::end() {
        return true;
    }

    bool contactAngleEvaluationTest::write() {
        if (time_.writeTime()) {
            contactAngles_.write();
            rf_.write();
            wisp_.write();
            rfMax_ = gMax(rf_);
            contactAnglemax_ = gMax (contactAngles_);
            reduce(contactAnglemin_, minOp<scalar>());
            if (Pstream::master())
            {
            // Create the output path directory
            fileName outputDir = "postProcessing";
            // Createe the directory
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr < OFstream > outputFilePtrMin;
            autoPtr < OFstream > outputFilePtrMax;
            autoPtr < OFstream > outputFilePtrRfMax;
            // Open the file in the newly created directory
            // TODO:  Maximum and min contact angle value -> might give erranious values if wisps are still there and detected as an interface cell.
            outputFilePtrMin.reset(new OFstream(outputDir/"contactangleMin"+".csv",
                                IOstreamOption(),
                                true));
            outputFilePtrMax.reset(new OFstream(outputDir/"contactangleMax"+".csv",
                                IOstreamOption(),
                                true));
            outputFilePtrRfMax.reset(new OFstream(outputDir/"rfMax"+".csv",
                                IOstreamOption(),
                                true));

            // Write stuff
            outputFilePtrMin() << time_.value()<< "," << contactAnglemin_<< endl;
            outputFilePtrMax() << time_.value()<< "," << contactAnglemax_<< endl;
            outputFilePtrRfMax() << time_.value()<< "," << rfMax_<< endl;
            //outputFilePtr() << time_.value();
        }

      }
      return true;
    }

}
// ************************************************************************* //
                                                 
