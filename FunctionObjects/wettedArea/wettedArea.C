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

#include "wettedArea.H"

#include "addToRunTimeSelectionTable.H"

#include <iostream>

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {

    defineTypeNameAndDebug(wettedArea, 0);
    addToRunTimeSelectionTable(functionObject, wettedArea, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    wettedArea::wettedArea(
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
        wettedArea_(0) {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool wettedArea::read(const dictionary & dict) {
        notImplemented("wettedArea::read(const dictionary&)");
        return execute();
    }

    bool wettedArea::execute() {
        calculateWettedArea();
        report();
        return true;
    }

    /// Estimates the wetted area (submerged area of a specific patch in the transported fluid)
    void wettedArea::calculateWettedArea() {
        wettedArea_ = 0;
        forAll(mesh_.boundary(), patch) {
            // Check if the patch is an alphaContactAngleTwoPhaseFvPatchScalarField
            const volScalarField::Boundary& abf = alpha1_.boundaryField();
            if (isA<alphaContactAngleTwoPhaseFvPatchScalarField>(abf[patch])) {
                // Get the patch name and ID number
                const word& patchName = mesh_.boundary()[patch].name();
                label patchID = mesh_.boundaryMesh().findPatchID(patchName);

                // Lookup the desired alpha values on the patch
                const fvPatchField<scalar>& alphaPatchField = alpha1_.boundaryField()[patchID];

                // Get the boundary mesh and the specific patch
                const fvBoundaryMesh& boundary = mesh_.boundary();
                const fvPatch& patch = boundary[patchID];

                // Get the area normal vectors of the patch
                const vectorField& patchSf = patch.Sf();
                // Calculate the wetted area by summing up the product of the area normal vectors and the alpha values
                forAll(patchSf, cellI) {
                    if (alphaPatchField[cellI] !=0) wettedArea_ += mag(patchSf[cellI]) * alphaPatchField[cellI];
                }
            }
	}
	// Synchronize the sum across processors
        reduce(wettedArea_, sumOp<scalar>());
    }



    //- Report at the console output / log file
    void wettedArea::report() {
        Info << endl << "Phase " << fieldName_ << " has " << (wettedArea_ * 1000000) <<
            " m^2 wetted area at time step" << time_.value() << endl;
    }

    bool wettedArea::start() {
        //execute();
        return true;
    }

    bool wettedArea::end() {
        return true;
    }

    bool wettedArea::write() {
    // Only the master processor writes to the output file
    if (Pstream::master()) {
        // Check if the current time step should be written to file
        if (time_.writeTime()) {
            // Create the output directory
            fileName outputDir = "postProcessing";
            mkDir(outputDir); // Create the directory

            // Open the output file in the newly created directory
            autoPtr<OFstream> outputFilePtr(new OFstream(outputDir / "wettedArea.csv", IOstreamOption(), true));

            // Write the current time and the wetted area (in square meters) to the file
            *outputFilePtr << time_.value() << "," << wettedArea_ * 1000000 << endl;
        }
    }
    return true;
}


}
// ************************************************************************* //