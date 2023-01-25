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

#include "contactAngle3D.H"

#include "addToRunTimeSelectionTable.H"

#include <iostream>

#include "vectorList.H"

#include <cmath>

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {

    defineTypeNameAndDebug(contactAngle3D, 0);
    addToRunTimeSelectionTable(functionObject, contactAngle3D, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    contactAngle3D::contactAngle3D(
            const word & name,
                const Time & time,
                    const dictionary & dict
        ):
        functionObject(name),
        time_(time),
        mesh_(time.lookupObject < fvMesh > (polyMesh::defaultRegion)),
        fieldName_(dict.lookup("phaseIndicator")),
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
        contactAngles_(
            IOobject(
                "contactAngles",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(
                "zero",
                dimless,
                0.0
            )
        ),
        wisp_(
            IOobject(
                "wisp",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(
                "zero",
                dimless,
                0.0
            )
        ),
        advPos_(
            IOobject(
                "advPos",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(
                "zero",
                dimless,
                0.0
            )
        ),
        recPos_(
            IOobject(
                "recPos",
                time.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar(
                "zero",
                dimless,
                0.0
            )
        ),
        contactAnglemin_(1000),
        contactAnglemax_(0),
        rfMax_(0) {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool contactAngle3D::read(const dictionary & dict) {
        notImplemented("contactAngle3D::read(const dictionary&)");
        return execute();
    }

    bool contactAngle3D::execute() {
        checkWisps();
        calContactAngle();
        report();
        return true;
    }

    //Calculates the contact angle in the interface cells, discarding the wisps
    void contactAngle3D::calContactAngle() {
        //Mesh connectivity
        const pointField & points = mesh_.points(); // Node coordinates
        const faceList & faces = mesh_.faces(); // Face to node
        const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion
        const volVectorField & C = mesh_.C(); // Cell center coordinates
        const auto & faceOwner = mesh_.faceOwner();
        contactAngles_ *= scalar(0);
        advPos_ *= scalar(0);
        recPos_ *= scalar(0);
        recPos_ += scalar(1000);
        contactAnglemax_ = 0;
        contactAnglemin_ = 1000;
        const volScalarField::Boundary & abf = alpha1_.boundaryField();

        scalar vAdv = 0.0; // velocity at advancing end
        scalar vRec = 0.0; // velocity at receding cell
        vector posAdv = vector::zero; // advnacing cell position
        vector posRec = vector::zero; // receding cell position
        scalar posAdvX = 0.0; // advnacing cell position x-component
        scalar posRecX = 0.0; // receding cell position x-component
       

        forAll(mesh_.boundary(), patch) {

            const word & patchName = mesh_.boundary()[patch].name(); // Boundary patch name

            if (isA < alphaContactAngleTwoPhaseFvPatchScalarField > (abf[patch])) {

                label patchID = mesh_.boundaryMesh().findPatchID(patchName);
                const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
                const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();
                int count = 0;
                int rcount = 0;
                rfMax_ = 0;
                const volVectorField & U = alpha1_.db().objectRegistry::lookupObject < volVectorField > ("U");
                const vectorField & Uw = U.boundaryField()[patchID];
                //loop over all faces of the boundary patch 

                forAll(mesh_.boundary()[patch], facei) {
                    // Get normals to the patch cell face
                    const vectorField & nHatf = mesh_.boundary()[patch].nf();
                    const label & face = boundaryMesh[patch].start() + facei; // Face index
                    int bCell = faceOwner[facei + mesh_.boundary()[patch].start()]; //cell ID

                    // check if the cell is an interface cell
                    if (alpha1_[bCell] < (1.0 - wispTol_) && alpha1_[bCell] > wispTol_) {
                        std::vector < double > signedDistancesOfPatchFace; // to store sgined distances of the interface path face
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
                        forAll(signedDistancesOfPatchFace, sd) {
                            auto ctr = (sd + 1) % signedDistancesOfPatchFace.size();
                            if (std::signbit(signedDistancesOfPatchFace[sd]) != std::signbit(signedDistancesOfPatchFace[ctr])) {
                                countSign++;
                            }
                        }
                        if (countSign != 0) // also counts for vertex intersection (2 for edge-edge intersection 4 for vertex-vertex, 2 for egde intersection), //3 for vertex-edge
                        {
                            isContactLine = true;
                        }
                        if (isContactLine) // detects the contact line
                        {
                            if (!wisp_[bCell]) // also checks if not a wisp 
                            {
                                if ((mag(plicNormalPatchInternalField[facei]) * mag(nHatf[facei])) != 0) {
                                    contactAngle = acos((plicNormalPatchInternalField[facei]) & nHatf[facei] / (mag(plicNormalPatchInternalField[facei]) * mag(nHatf[facei]))) * 180.0 / M_PI;
                                    contactAngles_[bCell] = contactAngle;
                                }

                                if (count == 0) {
                                    //Info << " I am in count 0 " << endl;
                                    posAdv = plicCentres_[bCell];
                                    
                                    posRec = plicCentres_[bCell];
                                    
                                    vAdv = mag(U[bCell]);
                                    vRec = mag(U[bCell]);
                                    posAdvX = posAdv.component(vector::X);
                                    posRecX = posRec.component(vector::X);
                                    //Info << posAdvX << " " << " " <<posRecX << " posAdv in cont 0" << endl;
                                    if (plicNormals_[bCell].component(vector::X) < 0) // this is hard coded for the sliding droplet which slides in the +x-direction
                                    //Can change it to the check for magnitude of the velocity
                                    {
                                        advPos_ *= scalar(0);
                                        advPos_[bCell] = posAdvX;
                                    }
                                    else if (plicNormals_[bCell].component(vector::X) > 0) 
                                    {
                                        recPos_ *= scalar(0);
                                        recPos_ += scalar(1000);
                                        recPos_[bCell] = posRecX;
                                    }

                                    contactAnglemin_ = contactAngle;
                                    count++;
                                } else if (count > 0) {
                                    if (plicNormals_[bCell].component(vector::X) < 0) { 
                                        posAdv = plicCentres_[bCell]; // can take plic centre 
                                       // Info << posAdvX << " posAdv in cont 1-0 " << bCell<< endl;
                                        if (posAdvX < posAdv.component(vector::X)){
                                            posAdvX = posAdv.component(vector::X); 
                                           // Info << posAdv << " posAdv in cont 1-0 " << bCell<< endl;
                                            vAdv = mag(U[bCell]);
                                            advPos_ *= scalar(0);
                                            advPos_[bCell] = posAdvX;
                                            //advPos_[bCell] = 1.0;
                                        }
                                    } else if (plicNormals_[bCell].component(vector::X) > 0) { 
                                        posRec = plicCentres_[bCell];
                                       // Info << posRec << " posAdv in cont 1-0 " << bCell<< endl;
                                        if (posRecX > posRec.component(vector::X)){
                                            posRecX = posRec.component(vector::X);
                                            //Info << posRecX << " posRec in cont 1-0 " << bCell << endl;
                                            vRec = mag(U[bCell]);
                                            recPos_ *= scalar(0);
                                            recPos_ += scalar(1000);
                                            recPos_[bCell] = posRecX;
                                            //recPos_[bCell] = -1.0;
                                        }

                                    }
                                }
                                if (contactAngle < contactAnglemin_) {

                                    contactAnglemin_ = contactAngle;
                                }
                            }
                        }
                    }
                }
            }
        }
        // reduce(vAdv, maxOp < scalar > ());
        //     reduce(vRec, minOp < scalar > ());
        //     Info << vAdv << " vadv " << vRec << " vRec " << endl;
    }

    //Check if an interface cell is a wisp. Marks all the wisp as 1 for visualization
    // void contactAngle3D::checkWisps() {
    //     wisp_ *= scalar(0);
    //     const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion 
    //     const auto & patches = mesh_.boundary();
    //     const auto & faceOwner = mesh_.faceOwner();
    //     const volScalarField::Boundary & abf = alpha1_.boundaryField();

    //     forAll(patches, patch) {
    //         const word & patchName = patches[patch].name(); // Boundary patch name
    //         if (isA < alphaContactAngleTwoPhaseFvPatchScalarField > (abf[patch])) {
    //             label patchID = boundaryMesh.findPatchID(patchName);
    //             const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
    //             const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();

    //             forAll(patches[patch], facei) {
    //                 const label & ownerCellId = faceOwner[facei + patches[patch].start()]; //cell ID
    //                 //if (plicNormalPatchInternalField[facei] != vector::zero)
    //                 if (alpha1_[ownerCellId] < (1.0 - alphaTol_) && alpha1_[ownerCellId] > alphaTol_) {
    //                     const labelList & faces = boundaryMesh[patch].faceFaces()[facei];
    //                     bool isAWisp = true;
    //                     forAll(faces, face) {
    //                         const auto & faceNId = faceOwner[faces[face] + patches[patch].start()]; //cell ID
    //                         if (alpha1_[faceNId] < (1.0 - wispTol_)  && alpha1_[faceNId] > wispTol_) { //full cell
    //                             isAWisp = false;
    //                             break;
    //                         }
    //                     }
    //                     if (isAWisp) {
    //                         wisp_[ownerCellId] = 1.0;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // Checks if an interface cell is a wisp and marks all the wisps as 1 for visualization
    void contactAngle3D::checkWisps() {
        wisp_ *= scalar(0); // Reset the wisp field
        const polyBoundaryMesh & boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion 
        const auto & patches = mesh_.boundary();
        const auto & faceOwner = mesh_.faceOwner();
        const volScalarField::Boundary & abf = alpha1_.boundaryField();

        forAll(patches, patch) {
            // Check if the patch is a boundary with contact angle boundary condition
            const word & patchName = patches[patch].name(); // Boundary patch name
            if (isA <alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patch])) {
                label patchID = boundaryMesh.findPatchID(patchName);
                const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchID].patchInternalField();
                const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchID].patchInternalField();

                forAll(patches[patch], facei) {
                    // Get the owner cell ID of the current face
                    const label & ownerCellId = faceOwner[facei + patches[patch].start()];

                    // Check if the owner cell is an interface cell
                    if (alpha1_[ownerCellId] < (1.0 - alphaTol_) && alpha1_[ownerCellId] > alphaTol_) {
                        const labelList & faces = boundaryMesh[patch].faceFaces()[facei];
                        bool isAWisp = true; // Assume that the cell is a wisp until proven otherwise
                        forAll(faces, face) {
                            // Get the owner cell ID of the face
                            const auto & faceNId = faceOwner[faces[face] + patches[patch].start()];
                            // Check if the cell is a non-wisp cell
                            if (alpha1_[faceNId] < (1.0 - wispTol_) && alpha1_[faceNId] > wispTol_) { //full cell
                                isAWisp = false;
                                break;
                            }
                        }
                        if (isAWisp) {
                            wisp_[ownerCellId] = 1.0; // Mark the cell as a wisp
                        }
                    }
                }
            }
        }
    }


    //- Report at the console output / log file
    void contactAngle3D::report() {
    }

    bool contactAngle3D::start() {
        //execute();
        return true;
    }

    bool contactAngle3D::end() {
        return true;
    }

    bool contactAngle3D::write() {
        if (time_.writeTime()) {
            contactAngles_.write();
            scalar advPosVal = gMax(advPos_);
            scalar recPosVal = gMin(recPos_);
            advPos_.write();
            recPos_.write();
            wisp_.write();
            // Create the output path directory
            fileName outputDir = "postProcessing";
            // Create the directory
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr < OFstream > outputFilePtrCA;
            autoPtr < OFstream > outputFilePtrCLPos;
            // Open the file in the newly created directory
            outputFilePtrCA.reset(new OFstream(outputDir / "contactAngle" + ".csv",
                IOstreamOption(),
                true));
            outputFilePtrCLPos.reset(new OFstream(outputDir / "contactLinePos" + ".csv",
                IOstreamOption(),
                true));

            // Get the maximum and minimum contact angles
            contactAnglemax_ = gMax(contactAngles_);
            reduce(contactAnglemin_, minOp < scalar > ());
            

            if (Pstream::master()) {
                // Write the contact angle data to the file
                outputFilePtrCA() << time_.value() << "," << contactAnglemin_ << "," << contactAnglemax_<< endl;
                // Write the advancing and receding contact line positions to the file
                outputFilePtrCLPos() << time_.value() << "," << recPosVal<<  "," <<advPosVal << endl;
            }
        }
        return true;
    }

}


    // bool contactAngle3D::write() {
    //     if (time_.writeTime()) {
    //         //reduce(advPos_, maxOp < scalar > ());
    //         contactAngles_.write();
    //         scalar advPosVal = gMax(advPos_);
    //         scalar recPosVal = gMin(recPos_);
    //         advPos_.write();
    //         recPos_.write();
    //         wisp_.write();
    //         fileName outputDir = "postProcessing";
    //         mkDir(outputDir);
    //         // File pointer to direct the output to
    //         autoPtr < OFstream > outputFilePtrCA;
    //         //autoPtr < OFstream > outputFilePtrMax;
    //         //autoPtr < OFstream > outputFilePtrRfMax;
    //         autoPtr < OFstream > outputFilePtrCLPos;
    //         outputFilePtrCA.reset(new OFstream(outputDir / "contactAngle" + ".csv",
    //             IOstreamOption(),
    //             true));
    //         // outputFilePtrMax.reset(new OFstream(outputDir / "contactangleMax" + ".csv",
    //         //     IOstreamOption(),
    //         //     true));
    //         // outputFilePtrRfMax.reset(new OFstream(outputDir / "rfMax" + ".csv",
    //         //     IOstreamOption(),
    //         //     true));
    //         outputFilePtrCLPos.reset(new OFstream(outputDir / "contactLinePos" + ".csv",
    //             IOstreamOption(),
    //             true));

    //         contactAnglemax_ = gMax(contactAngles_);

    //         reduce(contactAnglemin_, minOp < scalar > ());

    //         if (Pstream::master()) {
    //             // Open the file in the newly created directory
    //             // TODO:  Maximum and min contact angle value -> might give erranious values if wisps are still there and detected as an interface cell.
    //             // Write stuff
    //             outputFilePtrCA() << time_.value() << "," << contactAnglemin_ << "," << contactAnglemax_<< endl;
    //             //outputFilePtrMax() << time_.value() << "," << contactAnglemax_ << endl;
    //             //outputFilePtrRfMax() << time_.value() << "," << rfMax_ << endl;
    //             outputFilePtrCLPos() << time_.value() << "," << recPosVal<<  "," <<advPosVal << endl;
            
// ************************************************************************* //