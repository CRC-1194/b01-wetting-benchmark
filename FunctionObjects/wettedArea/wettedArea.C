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
        patchName_(dict.lookup("patchName")),
        alpha1_(
            mesh_.lookupObject < volScalarField >
            (
                fieldName_
            )
        ),
        wettedArea_(0)
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool wettedArea::read(const dictionary & dict) {
        notImplemented("wettedArea::read(const dictionary&)");
        return execute();
    }

    bool wettedArea::execute() {
        calcWettedArea();
        report();
        return true;
    }

    //Estimates the wetted area (submerged area of a specific patch in the trasported fluid)
    void wettedArea::calcWettedArea() {
        wettedArea_ = 0;
        //get the patch ID number
        label patchID = mesh_.boundaryMesh().findPatchID(patchName_);

        //Lookup the desired alpha values on the patch
        const fvPatchField < scalar > & alphaPatchField = alpha1_.boundaryField()[patchID];

        //- Get boundary mesh
        const fvBoundaryMesh & boundary = mesh_.boundary();
        //- Get desired path mesh from boundary mesh
        const fvPatch & patch = boundary[patchID];
        //- Get Area normal vectors of the desired patch
        const vectorField & patchSf = patch.Sf();

        //- Calculation of the wetted area
        forAll(patchSf, cellI) {
            wettedArea_ += mag(patchSf[cellI]) * alphaPatchField[cellI];
        }

        // Sync sum across processors
        reduce(wettedArea_, sumOp<scalar>());
    }

    //- Report at the console output / log file
    void wettedArea::report() {
        Info << endl << "Phase " << fieldName_ << " has " << (wettedArea_*1000000) <<
            " m^2 wetted area at time step" << time_.value()<< endl;
    }

    bool wettedArea::start() {
        //execute();
        return true;
    }

    bool wettedArea::end() {
        return true;
    }

    bool wettedArea::write() {
      if (Pstream::master())
      {
        if (time_.writeTime()) {
            // Create the output path directory
            fileName outputDir = "postProcessing";
            // Createe the directory
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr < OFstream > outputFilePtr;
            // Open the file in the newly created directory
            outputFilePtr.reset(new OFstream(outputDir/"wettedArea.csv",
                                IOstreamOption(),
                                true));

            // Write stuff
            outputFilePtr() << time_.value() << ","
                            << wettedArea_*1000000 << endl;
        }
      }
      return true;
    }

}
// ************************************************************************* //
