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

#include "contactAngleEvaluation.H"

#include "addToRunTimeSelectionTable.H"

#include <iostream>
#include "fvc.H"
#include "vectorList.H"

#include <cmath>

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"

#include "volMesh.H"
#include "volFields.H"

#include "mathematicalConstants.H"

#include <functional>

#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {

    defineTypeNameAndDebug(contactAngleEvaluation, 0);
    addToRunTimeSelectionTable(functionObject, contactAngleEvaluation, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    contactAngleEvaluation::contactAngleEvaluation(
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


    bool contactAngleEvaluation::read(const dictionary & dict) {
        notImplemented("contactAngleEvaluation::read(const dictionary&)");
        return execute();
    }

    bool contactAngleEvaluation::execute() {
        checkWisps();
        calContactAngle();
        report();
        return true;
    }

    // It calculates the maximum contact radius from all the contact points in 3D domain
    void contactAngleEvaluation::maxRf (const volScalarField& rf_)
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
    void contactAngleEvaluation::calRf(label bCell, const volVectorField& C, int &rcount)
    {
        // provides the distance between the central axis and the contact line cell centre
        label cellID = bCell;
        point origin = origin_;
        point contactLinePoint = C[cellID];
        point positionOnCentralAxis = {origin[0], origin[1], contactLinePoint[2]};
       
        if (mag((contactLinePoint - positionOnCentralAxis)) != 0){
            rf_[cellID] = (contactLinePoint - positionOnCentralAxis) & ((contactLinePoint - positionOnCentralAxis) / mag((contactLinePoint - positionOnCentralAxis)));
        }
    }


    //Calculates the contact angle in the interface cells, discarding the wisps
    void contactAngleEvaluation::calContactAngle() {
        //Mesh connectivity
        const volVectorField& C = mesh_.C();         // Cell center coordinates
        const auto& faceOwner = mesh_.faceOwner();
        contactAngles_*= scalar(0);
        rf_*= scalar(0);
	    const volScalarField::Boundary & abf = alpha1_.boundaryField();
        const fvBoundaryMesh& boundary = mesh_.boundary();


        forAll(boundary, patch)
        {
	        if (isA < alphaContactAngleTwoPhaseFvPatchScalarField > (abf[patch])) 
            {
                const word& patchName = boundary[patch].name();            // Boundary patch name
                const label patchIndex = boundary[patch].index();

                //patch face normal vectors
                const vectorField nf(boundary[patch].nf());
                
                //label patchID = mesh_.boundaryMesh().findPatchID(patchName);
                // Get patch internal fields of normals and centers
                const vectorField plicNormalPatchInternalField = plicNormals_.boundaryField()[patchIndex].patchInternalField();
                const vectorField plicCentrePatchInternalField = plicCentres_.boundaryField()[patchIndex].patchInternalField();
                
                // surfaceVectorField to compute theta
                surfaceVectorField plicNormalsf(fvc::interpolate(plicNormals_));
                forAll(plicNormalsf.boundaryFieldRef()[patch],i)
                {
                    const label celli = boundary[patch].faceCells()[i];
                    vector n = plicNormals_[celli];
                    if(mag(n) != 0)
                    {
                        n /= mag(n);
                        plicNormalsf.boundaryFieldRef()[patch][i] = n;
                    }
                }

                fvsPatchVectorField& plicNormalPatchField = plicNormalsf.boundaryFieldRef()[patch];
                
                // Compute the contact angles at the wall.
                tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(plicNormalPatchField & nf));
                scalarField& thetaf = thetafTmp.ref();

                int rcount = 0;
                rfMax_ = 0;
                
                //loop over all faces of the boundary patch 
                forAll(thetaf, faceI)
                {
                        double contactAngle = 0.0;
                        if(hasContactLine(faceI, boundary[patch]))// detects the contact line
                        { 
                            int cellI = faceOwner[faceI + mesh_.boundary()[patch].start()]; // global cell index
                            if(!wisp_[cellI]) // also checks if not a wisp  
                            {
                                calRf(cellI, C, rcount);
                                if((mag (plicNormals_[cellI]) * mag(nf[faceI])) !=0)
                                {
                                    contactAngle = acos((plicNormals_[cellI]) & nf[faceI] / (mag (plicNormals_[cellI]) * mag(nf[faceI])) ) * 180.0 / M_PI;
                                    contactAngles_[cellI] = contactAngle;
				                }
                            }
                        }
                }

            }
        }
    }

//Check if an interface cell is a wisp. Marks all the wisp as 1 for visualization
    void contactAngleEvaluation::checkWisps(){
        wisp_*= scalar(0);
        const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh(); //boundary mesh informataion 
        const auto& patches = mesh_.boundary();
        const auto& faceOwner = mesh_.faceOwner();
	    const volScalarField::Boundary & abf = alpha1_.boundaryField();

        forAll(patches, patch)
        {
            const word& patchName = patches[patch].name();            // Boundary patch name
            if (isA < alphaContactAngleTwoPhaseFvPatchScalarField > (abf[patch])) 
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

    // Check if the patch face has a contact line: based on signed distance calculations
    bool contactAngleEvaluation::hasContactLine(label faceI, const fvPatch& patch) const
    {
        // Look up PLIC normals and positions. 
        const auto& db = alpha1_.db(); 

        const auto normalsName = IOobject::groupName
        (
            "interfaceNormal", 
            alpha1_.internalField().group()
        );
        const auto centresName = IOobject::groupName
        (
            "interfaceCentre", 
            alpha1_.internalField().group()
        );

        bool hasNormals = db.found(normalsName);
        if (!hasNormals)
        {
            // This BC is updated before interface reconstruction.
            // Do nothing if PLIC fields are not available in the registry. 
            return false;
        }

        bool hasCentres = db.found(centresName);
        if (!hasCentres)
        {
            // This BC is updated before interface reconstruction.
            // Do nothing if PLIC fields are not available in the registry. 
            return false;
        }

        const volVectorField& interfaceNormal = 
            db.lookupObject<volVectorField>(normalsName);

        const volVectorField& interfaceCentre = 
            db.lookupObject<volVectorField>(centresName);


        // Get patch fields for interface normals and centers
        const label patchIndex = patch.index();
        const auto& pInterfaceNormals = interfaceNormal.boundaryField()[patchIndex];
        const auto& pInterfaceCentres = interfaceCentre.boundaryField()[patchIndex];

        // Get patch internal fields of normals and centers
        const auto pInternalNormalsTmp = pInterfaceNormals.patchInternalField();
        const auto& pInternalNormals = pInternalNormalsTmp.cref(); 
        const auto pInternalCentresTmp = pInterfaceCentres.patchInternalField();
        const auto& pInternalCentres = pInternalCentresTmp.cref(); 

        const vector& cellInterfaceNormal = pInternalNormals[faceI];
        const vector& cellInterfaceCentre = pInternalCentres[faceI];

        const auto& mesh = interfaceNormal.mesh();
        const auto& meshPoints = mesh.points();
        const auto& meshFaces = mesh.faces();
        const auto& thisFace = meshFaces[patch.start() + faceI];

        // Get face points. 
        for(auto pointI = 0; pointI < (thisFace.size() - 1); ++pointI)
        {
            // Compute the signed distance of the first point.
            const point& firstFacePoint = meshPoints[thisFace[pointI]];
            const scalar firstDist = (firstFacePoint - cellInterfaceCentre) & 
                cellInterfaceNormal;

            // Compute the signed distance of the second point.
            const point& secondFacePoint = meshPoints[thisFace[pointI + 1]];
            const scalar secondDist = (secondFacePoint - cellInterfaceCentre) & 
                cellInterfaceNormal;

            if (firstDist * secondDist < 0)
            {
                return true;
            }
        }

        return false;
    }

    //- Report at the console output / log file
    void contactAngleEvaluation::report() {}

    bool contactAngleEvaluation::start() {
        return true;
    }

    bool contactAngleEvaluation::end() {
        return true;
    }

    bool contactAngleEvaluation::write() {
        if (time_.writeTime()) {
            contactAngles_.write();
            rf_.write();

            //wisp_.write();
            rfMax_ = gMax(rf_);
            //contactAnglemax_ = gMax (contactAngles_);
            //reduce(contactAnglemin_, minOp<scalar>());
            if (Pstream::master())
            {
              //  std::ofstream outputFileMin("postProcessing/contactangleMin.csv", std::ios::app);
              //  std::ofstream outputFileMax("postProcessing/contactangleMax.csv", std::ios::app);
                std::ofstream outputFileRfMax("postProcessing/rfMax.csv", std::ios::app);
                // TODO:  Maximum and min contact angle value -> might give erranious values if wisps are still there and detected as an interface cell.

                // Write stuff
               // outputFileMin << time_.value()<< "," << contactAnglemin_<< nl;
               // outputFileMax << time_.value()<< "," << contactAnglemax_<< nl;
                outputFileRfMax << time_.value()<< "," << rfMax_<< nl;
            
            }
        }
        return true;
    }

}
// ************************************************************************* //
                                                 
