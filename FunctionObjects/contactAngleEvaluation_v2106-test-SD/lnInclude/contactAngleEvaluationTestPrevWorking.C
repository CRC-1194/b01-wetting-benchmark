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
    void contactAngleEvaluationTest::calRf(label bCell, const volVectorField& C, int &rcount)
    {
        // provides the distance between the central axis and the contact line cell centre
        label cellID = bCell;
        point origin = origin_;
        point contactLinePoint = C[cellID ];
        point positionOnCentralAxis = {origin[0], origin[1], contactLinePoint[2]};
        //label celli = mesh_.findNearestCell(positionOnCentralAxis);
        //point pOnCentralAxis = C[celli];
        if (mag((contactLinePoint - positionOnCentralAxis)) != 0){
            rf_[cellID] = (contactLinePoint - positionOnCentralAxis) & ((contactLinePoint - positionOnCentralAxis) / mag((contactLinePoint - positionOnCentralAxis)));
           // if (Pstream::master())
            //std::cout << rf_[cellID] << " rf and count " << rcount << " cell id " << cellID << std::endl;
            //std::cout << mag((contactLinePoint - pOnCentralAxis)) << " mag and cell ID " << cellID << std::endl; 
        }
        //if(rcount==0)
        //{
          //  rfMax_ +=  rf_[cellID];
            //rcount++;
        //}
        //else
       // {
         //   if(rf_[cellID] > rfMax_) {
           //    rfMax_ +=  rf_[cellID];
               //rcount++;
           // } 
           // rcount++;   
       // } 
    }



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

            if (patchName == "bottom")
            {
                label patchID = mesh_.boundaryMesh().findPatchID(patchName);
                const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
                const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();
                int count=0;
                int rcount = 0;
                rfMax_ = 0;
                //std::cout << mesh_.boundary()[patch].size() << " size " << std::endl;
                 //loop over all faces of the boundary patch 
                forAll(mesh_.boundary()[patch], facei)
                {
                    // Get normals to the patch cell face
                    const vectorField& nHatf = mesh_.boundary()[patch].nf();           
                    //const label& bCell = boundaryMesh[patch].faceCells()[facei];    // Boundary cell index / Owner cell index
                    const label& face = boundaryMesh[patch].start() + facei;        // Face index
                    int bCell = faceOwner[facei + mesh_.boundary()[patch].start()]; //cell ID
                    
                    //if the patch cell has an interface
                    //if (plicNormalPatchInternalField[facei] != vector::zero)
                    if(alpha1_[bCell] < (1.0-wispTol_) && alpha1_[bCell] > wispTol_ )
                    {

                        // if(alpha1_[bCell] > alphaTol_)
                        // {
                        //     if(alpha1_[bCell] < (1.0-alphaTol_))
                        //     {
                                //if (Pstream::master())
                                //std::cout << alpha1_[bCell] << " alpha value " << std::endl;
                                std::vector <double> signedDistancesOfPatchFace; // to store sgined distances of the interface path face
                                bool isContactLine = false;
                                int countSign = 0;
                                double contactAngle = 0.0;
                                

                                forAll(faces[face], vertexID) //patch face
                                {
                                    point vertex = points[faces[face][vertexID]]; // vertices of a patch face with interface
                                    double sd_ = plicNormalPatchInternalField[facei] & (plicCentrePatchInternalField[facei] - vertex);
                                    signedDistancesOfPatchFace.push_back(sd_);
                                    //if (Pstream::master())
                                   // std::cout << " PLIC components: Centre " <<  plicCentrePatchInternalField[facei][0] << " "<<  plicCentrePatchInternalField[facei][1] << " "<<  plicCentrePatchInternalField[facei][2] << " "
                                     //           <<  plicNormalPatchInternalField[facei][0] << " " << plicNormalPatchInternalField[facei][1] << " " << plicNormalPatchInternalField[facei][2] << std::endl;
                                    //if (Pstream::master())
                                    //std::cout << " Cell ID "  << bCell << " Point " << vertex[0] << " " << vertex[1] << " " << vertex[2] << " SD: " << sd_ <<std::endl;
                                    
                                }    
                                forAll(signedDistancesOfPatchFace, sd)
                                {
                                    auto ctr = (sd + 1) % signedDistancesOfPatchFace.size();
                                    if (std::signbit(signedDistancesOfPatchFace[sd]) != std::signbit(signedDistancesOfPatchFace[ctr]))
                                    {
                                        countSign++;
                                    }
                                    //if (Pstream::master())
                                    //std::cout << " inside sd calculation " << bCell<< std::endl;
                                }
                                if(countSign!=0) // also counts for vertex intersection (2 for edge-edge intersection 4 for vertex-vertex, 2 for egde intersection), //3 for vertex-edge
                                {
                                    isContactLine = true;
                                    calRf(bCell, C, rcount);
                                    //if (Pstream::master())
                                    //std::cout << " inside count sign calculation " << bCell<< std::endl;
                                }
                                if(isContactLine){
                                    if(!wisp_[bCell])
                                    {
                                        //if (Pstream::master())
                                        //std::cout << " inside contact angle calculation " << bCell<< std::endl;
                                        contactAngle = acos((plicNormalPatchInternalField[facei]) & nHatf[facei] / (mag (plicNormalPatchInternalField[facei]) * mag(nHatf[facei])) ) * 180.0 / M_PI;
                                        contactAngles_[bCell] = contactAngle;
                                        //std::cout << " CA: cell ID " << bCell << " with plic normal " << plicNormalPatchInternalField[facei][0] << " ";
                                        //std::cout << plicNormalPatchInternalField[facei][1] << " " << plicNormalPatchInternalField[facei][2] << " and alpha " << alpha1_[bCell] << " cont A " << contactAngles_[bCell] << std::endl;
                                        
//                                        const labelList& faces = boundaryMesh[patch].faceFaces()[facei];
                                        //forAll(faces, face)
                                       // {
                                        //    const auto& faceNId = faceOwner[faces[face] + mesh_.boundary()[patch].start()]; //cell ID
                                            //if(isAWisp){
                                          //      std::cout << " CA: cell ID " << bCell << " with nighbour cells " << faceNId << " ";
                                            //    std::cout << " CA: alpha " <<alpha1_[faceNId] << " cont A " << contactAngles_[bCell] << std::endl;
                                            //}
                                            // if(alpha1_[faceNId] < (1.0- alphaTol_) ){ //full cell
                                            //     isAWisp = false; 
                                            //     std::cout << " cell ID " << ownerCellId << " with nighbour cells " << faceNId << " ";
                                            //     std::cout << " alpha " <<alpha1_[faceNId] << " cont A " << contactAngles_[ownerCellId] << std::endl;
                                            //     break;
                                            // }
                                       // }

                                        if(count==0)
                                        {
                                            contactAnglemax_ = contactAngle;
                                            contactAnglemin_ = contactAngle;
                                            count++;
                                        }
                                        if(contactAngle > contactAnglemax_)
                                        {
                                            contactAnglemax_ = contactAngle;
                                        }
                                        if(contactAngle < contactAnglemin_)
                                        {
                                            contactAnglemin_ = contactAngle;
                                        }
                                        //if (Pstream::master())
                                        //std::cout << contactAngle << " contact angle of the cell ID : " << bCell << std::endl;
                                    }
                                }
                                //rfMax_ /= (rcount+1);
                        //     }
                        // }
                    }
                }
                // Sync sum across processors
            //reduce(rfMax_, sumOp<scalar>());
            }
        }
    }

    void contactAngleEvaluationTest::checkWisps(){
       wisp_*= scalar(0);
        const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion 
        const auto& patches = mesh_.boundary();
        const auto& faceOwner = mesh_.faceOwner();

        forAll(patches, patch)
        {
            const word& patchName = patches[patch].name();            // Boundary patch name
            if(patchName == "bottom")
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
                            if(isAWisp){
                                //std::cout << " Wisp: cell ID " << ownerCellId << " with nighbour cells " << faceNId << " ";
                                //std::cout << " Wisp: alpha " <<alpha1_[faceNId] << " cont A " << contactAngles_[ownerCellId] << std::endl;
                            }
                            if(alpha1_[faceNId] < (1.0- wispTol_) ){ //full cell
                                isAWisp = false; 
                                ///std::cout << " cell ID " << ownerCellId << " with nighbour cells " << faceNId << " ";
                                //std::cout << " alpha " <<alpha1_[faceNId] << " cont A " << contactAngles_[ownerCellId] << std::endl;
                                break;
                            }
                        }
                        if(isAWisp){
                            //std::cout << " WISP cll ID " << plicNormalPatchInternalField[ownerCellId][0] << " " << plicNormalPatchInternalField[ownerCellId][1] << " " << plicNormalPatchInternalField[ownerCellId][2]<<std::endl;
                            //std::cout << " Wisp: cell ID " << ownerCellId << " with contact angle  "  << contactAngles_[ownerCellId] << " alpha value " << alpha1_[ownerCellId] << std::endl;
                            wisp_[ownerCellId] =1.0;
                        }
                        //std::cout << std::endl;
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
        //calContactAngle(); 
      //if (Pstream::master())
      //{
        if (time_.writeTime()) {
            //std::cout << " ->  " << contactAnglemax_  << " " << contactAnglemin_<< std::endl;
            contactAngles_.write();
            rf_.write();
            wisp_.write();
            maxRf(rf_);   
            
            //Todo export to csv
            // Create the output path directory
            fileName outputDir = "postProcessing";
            // Createe the directory
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr < OFstream > outputFilePtrMin;
            autoPtr < OFstream > outputFilePtrMax;
            autoPtr < OFstream > outputFilePtrRfMax;
            // Open the file in the newly created directory
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
            /*forAll(contactAngle, ca)
                outputFilePtr() << ","
                            << contactAngle[ca];
            outputFilePtr() << endl;
            */
           
        }
        
      //}
     
      return true;
    }

}
// ************************************************************************* //
