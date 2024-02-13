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
    employs the empirical Kistler model for dynamic contact angle 
    ((taken from Dirk Drunding and Anja Lippert dissertation)) i.e.,
    
    theta_d = f_Hoff (Ca + f_Hoff^(-1) (theta_e));
    f_Hoff(x) = arccos(1-2*tanh(5.16(Ca / (1 + 1.31*Ca^0.99))^0.706))

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt
\*---------------------------------------------------------------------------*/

#include "kistler.H"

#include "addToRunTimeSelectionTable.H"

#include "fvPatchFieldMapper.H"

#include "volMesh.H"

#include "volFields.H"

#include "mathematicalConstants.H"

#include <cmath>

#include <functional>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kistler::
kistler
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0)
{}


Foam::kistler::
kistler
(
    const kistler& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_)
{}


Foam::kistler::
kistler
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(dict.get<scalar>("theta0"))
{
    evaluate();
}


Foam::kistler::
kistler
(
    const kistler& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_)
{}


Foam::kistler::
kistler
(
    const kistler& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::kistler::theta
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
                        180 / constant::mathematical::pi * kistler_model(theta0_, Ca[faceI]),
                        scalar(180)
                        );
        }
    }

    return thetafTmp;
}

// Newton's method for root finding
double Foam::kistler::newton_method(const std::function<double(double)>& func, const std::function<double(double)>& deriv, scalar initial_guess) const
{
    scalar x = initial_guess;
    scalar h = func(x) / deriv(x);
    while (std::abs(h) >= 1e-6) {
        h = func(x) / deriv(x);
        x = x - h;
    }
    return x;
}

//Calculate the inverse Hoffmann function
double Foam::kistler::hoffman_function_inv(scalar theta0_) const
{
    scalar theta0_rad = theta0_* M_PI / 180.0;
    scalar right_side = std::pow((std::atanh((1 - std::cos(theta0_rad)) / 2) / 5.16), 1 / 0.706);

    auto equation = [right_side](double y) -> double {
        return y - (1.31 * right_side * std::pow(y, 0.99)) - right_side;
    };

    auto derivative = [right_side](double y) -> double {
        return 1 - (1.31 * 0.99 * right_side * std::pow(y, -0.01));
    };

    scalar initial_guess = 0.1;
    double reValue =newton_method(equation, derivative, initial_guess);
    return reValue;
}


// Kistler model
double Foam::kistler::kistler_model(scalar theta0_, scalar Ca)  const 
{
    double invHoff = hoffman_function_inv(theta0_);
    double x = invHoff + mag(Ca);
    double theta_d = std::acos(1 - 2 * std::tanh(5.16 * std::pow(x / (1 + 1.31*std::pow(x, 0.99)), 0.706)));
    return theta_d;
}

void Foam::kistler::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntry("theta0", theta0_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        kistler
    );
}


// ************************************************************************* //
