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

#include "updateFlux.H"

#include "addToRunTimeSelectionTable.H"

#include <iostream>

#include "vectorList.H"
#include <list>
#include <cmath>
#include "fvCFD.H"

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {

    defineTypeNameAndDebug(updateFlux, 0);
    addToRunTimeSelectionTable(functionObject, updateFlux, dictionary);

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    updateFlux::updateFlux(
            const word & name,
                const Time & time,
                    const dictionary & dict
        ):
        functionObject(name),
        time_(time),
        mesh_(time.lookupObject < fvMesh > (polyMesh::defaultRegion)),
        fieldName_(dict.lookup("phaseIndicator")),
        acceleration_(dict.lookup("acceleration")),
        initialHeight_(dict.lookup("initialHeight")),
        alpha1_(
            mesh_.lookupObject < volScalarField >
            (
                fieldName_
            )
        ),
        updatedCells_(
            IOobject(
                "updatedCells",
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
        updateFlux_(0) {
            calculateupdateFlux();
        }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    bool updateFlux::read(const dictionary & dict) {
        notImplemented("updateFlux::read(const dictionary&)");
        return execute();
    }

    bool updateFlux::execute() {
        calculateupdateFlux();
        report();
        return true;
    }

    /// Estimates the wetted area (submerged area of a specific patch in the transported fluid)
    void updateFlux::calculateupdateFlux() {

        // get the initial height -> done via reading the controlDict
        // get the acceleration of the plate -> done via reading the controlDict
        // get the current time t of the simulation
        // Update the velocity v_c = a*t
        // Update the next position h(t) -> corresponding to the height below which all cell field values (U, flux) are to be updated
        // h(t) = at^2 / 2.0

        //Current time
        scalar currentTime = time_.value();
        scalar deltaT = time_.deltaTValue();

        // Update the velocity v_c = a*t -> uniform for all cells below h(t)
        Foam::vector updatedVelocity = acceleration_ * currentTime;

        // Update the next position h(t) -> corresponding to the height below which all cell field values (U, flux) are to be updated
        // h(t) = at^2 / 2.0
        Foam::vector updatedHeight = initialHeight_ +  (acceleration_ * sqr(currentTime) / 2.0);

        // get all cell IDs that are below h(t) -> cellsToUpdate()
        updatedCells_*=scalar(0);

        //get the velocity field
        Foam::volVectorField& U  = const_cast<Foam::volVectorField&> (alpha1_.db().objectRegistry::lookupObject < volVectorField > ("U"));
  
        forAll(mesh_.C(), cellI)
        {
            if((alpha1_[cellI] > (1e-8)) && (mesh_.C()[cellI].component(vector::Z) <= updatedHeight.component(vector::Z)))
            {
                updatedCells_[cellI] = 1.0;
                U[cellI] = updatedVelocity;
            }
        }

        // Get the flux field
        Foam::surfaceScalarField& phi = const_cast<Foam::surfaceScalarField&>  (alpha1_.db().objectRegistry::lookupObject<Foam::surfaceScalarField>("phi"));

        //phi = fvc::flux(U);
        phi = fvc::interpolate(U) & mesh_.Sf();

        if (time_.writeTime()){
            U.write();
            phi.write();
        }
        U.correctBoundaryConditions();
    }



    //- Report at the console output / log file
    void updateFlux::report() {}

    bool updateFlux::start() {
        //execute();
        return true;
    }

    bool updateFlux::end() {
        return true;
    }

    bool updateFlux::write() {
        if (time_.writeTime()){
            updatedCells_.write();
        }
   
    return true;
    }

}
// ************************************************************************* //