/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    A dynamic alphaContactAngle scalar boundary condition
    employs the hydrodynamic coxVoinov model for dynamic contact angle 
    ((taken from Dirk Drunding and Anja Lippert dissertation)) i.e.,
    
    // thetad = (beta*  Ca + theta0^3)^(1/3) for \theta<135° (Dirk Grunding dissertation)
    // beta = ln(x/L) -> x is the macro length scale and L is the micro length scale

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt
\*---------------------------------------------------------------------------*/

#include "coxVoinov.H"

#include "addToRunTimeSelectionTable.H"

#include "fvPatchFieldMapper.H"

#include "volMesh.H"

#include "volFields.H"

#include "mathematicalConstants.H"

#include <cmath>

#include <functional>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coxVoinov::
coxVoinov
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    beta_(0.0)
{}


Foam::coxVoinov::
coxVoinov
(
    const coxVoinov& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    beta_(gcpsf.beta_) 
{}


Foam::coxVoinov::
coxVoinov
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(dict.get<scalar>("theta0")),
    beta_(readScalar(dict.lookup("beta"))) 
{
    evaluate();
}


Foam::coxVoinov::
coxVoinov
(
    const coxVoinov& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    beta_(gcpsf.beta_)
{}


Foam::coxVoinov::
coxVoinov
(
    const coxVoinov& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    beta_(gcpsf.beta_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::coxVoinov::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{

    const vectorField nf(patch().nf());

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);
    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculated the component of the velocity parallel to the wall
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;
    // Calculate Uwall resolved normal to the interface parallel to
    // the interface
    scalarField uwall(nWall & Uwall);

    //Calculaton of Capillary number
    const label patchi = this -> patch().index();

    const volScalarField & nu1 =
        this -> db().objectRegistry::lookupObject < volScalarField > ("nu1");

    const dictionary & transportProperties =
        this -> db().objectRegistry::lookupObject < IOdictionary >
        (
            "transportProperties"
        );

    word phase1Name(wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name(wordList(transportProperties.lookup("phases"))[1]);


    dimensionedScalar rho1(transportProperties.subDict(phase1Name).getScalar("rho"));
    dimensionedScalar sigmap(transportProperties.getScalar("sigma"));


    const fvPatchScalarField & nu1p = nu1.boundaryField()[patchi];
    
    //Capillary number fields
    scalarField Ca(nu1p * rho1.value() * uwall / sigmap.value());

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref();
    // For all boundary faces
    forAll(thetaf, faceI)
    {    
        // If we are in a contact-line cell
        if (mag(nHat[faceI]) > 0) 
        {
            thetaf[faceI] = min(
                                180 / constant::mathematical::pi * (pow(9.0*beta_ *  mag(Ca[faceI]) +
                                    pow(theta0_ * constant::mathematical::pi / 180, 3),
                                    0.3333333)),
                                scalar(180)
                            );
        }
    }

    return thetafTmp;
}

void Foam::coxVoinov::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coxVoinov
    );
}


// ************************************************************************* //
