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

#include "interpolationCellPoint.H"

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

// Newton's method for root finding
double Foam::kistler::newton_method(
    const std::function<double(double)>& func, 
    const std::function<double(double)>& deriv, 
    scalar initial_guess
) const
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
    scalar right_side = std::pow(
        (std::atanh((1 - std::cos(theta0_rad)) / 2) / 5.16), 1 / 0.706);

    auto equation = [right_side](double y) -> double {
        return y - (1.31 * right_side * std::pow(y, 0.99)) - right_side;
    };

    auto derivative = [right_side](double y) -> double {
        return 1 - (1.31 * 0.99 * right_side * std::pow(y, -0.01));
    };

    scalar initial_guess = 0.1;
    double reValue = newton_method(equation, derivative, initial_guess);
    return reValue;
}


// Kistler model
double Foam::kistler::kistler_model(
    scalar theta0_, 
    scalar Ca
)  const 
{
    double invHoff = hoffman_function_inv(theta0_);
    double x = invHoff + mag(Ca);
    double theta_d = std::acos(
        1 - 2 * std::tanh(5.16 * 
                          std::pow(x / (1 + 1.31*std::pow(x, 0.99)), 0.706)));
    return theta_d;
}

// Check if the patch face has a contact line: based on 
// signed distance calculations
bool Foam::kistler::hasContactLine(label faceI) const
{
    // Look up PLIC normals and positions. 
    const auto& db = this->db(); 

    const auto normalsName = 
        IOobject::groupName("interfaceNormal", this->internalField().group());
    const auto centresName = 
        IOobject::groupName("interfaceCentre", this->internalField().group());

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
    const fvPatch& patch = this->patch();
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kistler::kistler(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:   
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


Foam::kistler::kistler(
    const kistler& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


Foam::kistler::kistler
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    theta0_(dict.get<scalar>("theta0")),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{
    evaluate();
}


Foam::kistler::kistler
(
    const kistler& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


Foam::kistler::kistler
(
    const kistler& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    contactLineAngle_
    (
        IOobject
        (
            "clangle", 
            this->patch().boundaryMesh().mesh().time().timeName(),  
            this->patch().boundaryMesh().mesh(), 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        this->patch().boundaryMesh().mesh(), 
        dimensionedScalar("clangle", dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>Foam::kistler::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    //Patch Index
    const auto& patch = this->patch();
    const label patchIndex =patch.index();

    //patch face normal vectors
    const vectorField nf(patch.nf());
    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);
    // Normalise nWall
    nWall /= (mag(nWall) + SMALL);

    // Calculate contact line velocity relative to the 
    // wall velocity for moving walls. 
    vectorField Uwall(Up.patchInternalField() - Up);

    // Calculate component of the contact line velocity Uwal in the direction
    // of the interface normal tagential to the wall.
    scalarField uwall(nWall & Uwall);

    const dictionary& transportProperties =
    this->db().objectRegistry::lookupObject<IOdictionary>
    (
        "transportProperties"
    );

    // Choose the larger dynamic viscosity from available two phases.
    word phase1Name (wordList(transportProperties.lookup("phases"))[0]);
    word phase2Name (wordList(transportProperties.lookup("phases"))[1]);

    // Get constant phase-specific densities and kinematic viscosities.
    dimensionedScalar rho1(transportProperties.subDict(phase1Name).get<dimensionedScalar>("rho"));
    dimensionedScalar nu1c(transportProperties.subDict(phase1Name).get<dimensionedScalar>("nu"));

    dimensionedScalar rho2(transportProperties.subDict(phase2Name).get<dimensionedScalar>("rho"));
    dimensionedScalar nu2c(transportProperties.subDict(phase2Name).get<dimensionedScalar>("nu"));

    word nuName; 
    dimensionedScalar rho(rho1); 
    // If the dynamic viscosity of phase1 is larger
    if (rho1*nu1c > rho2*nu2c)
    {
        nuName = "nu1";
        rho = rho1;
    }
    else
    {
        nuName = "nu2";
        rho = rho2;
    }
    const volScalarField& nu =
        this->db().objectRegistry::lookupObject<volScalarField>(nuName);

    // Fetch the wall kinematic viscosity of phase1 
    const fvPatchScalarField& nup = nu.boundaryField()[patchIndex];
    scalarField muwall (nup*rho.value());

    // Fetch surface tension coefficient
    dimensionedScalar sigmap(transportProperties.get<dimensionedScalar>("sigma"));

    word fieldName = "alpha." + phase1Name;
    const volScalarField& alpha1_ =
        nu.mesh().lookupObject<volScalarField>(fieldName);

    const volScalarField::Boundary& abf = alpha1_.boundaryField();   
    // Lookup the desired alpha values on the patch
    const fvPatchField<scalar>& alphaPatchField = abf[patchIndex];
    
    // Capillary number will be updated later in the code 
    // based on (i) cell-centred velocity field
    // (ii) interpolated velocity field at interface-centres

    // Wall Capillary number
    scalarField Ca(muwall*uwall/sigmap.value());

    // 1. Get the PLIC centre and normals
    // Look up PLIC normals and positions. 
    const auto& db = this->db(); 

    const auto normalsName = 
        IOobject::groupName("interfaceNormal", this->internalField().group());
    const auto centresName = 
        IOobject::groupName("interfaceCentre", this->internalField().group());

    bool hasNormals = db.found(normalsName);
    bool hasCentres = db.found(centresName);
            
    // Cell-centred velocity field
    Foam::volVectorField& U  = const_cast<Foam::volVectorField&> (this -> db().objectRegistry::lookupObject < volVectorField > ("U"));

    // NOTE: Velocity field at the plic centre
    // is calculated using inverse distance interpolation scheme using U and cell and plic centres
    // This field is then used to calculate the Capillary number in this code
    // which is consistent with isoAdvection.C advection scheme
    
    // Uplic is used for advection: interpolated velocity at the PLIC centre
    Foam::volVectorField Uplic (U*scalar(0));

    // If there are no PLIC elements (i) calculate Ca  using cell-centred velocity field
    // Else calculate Ca using interpolated cell centred velocity to PLIC centre (using inverse distance weightings)
    if (!hasNormals || !hasCentres )
    {
        //Info << "not Interpolating the PLIC" << nl;
        // This BC is updated before interface reconstruction.

        // Calculated the component of the velocity parallel to the wall
        Uwall = Up.patchInternalField() - Up;
        Uwall -= (nf & Uwall) * nf;

        // Find the direction of the interface parallel to the wall
        nWall = nHat - (nf & nHat) * nf;

        // Normalise nWall
        nWall /= (mag(nWall) + SMALL);

        // Calculate Uwall resolved normal to the interface parallel to
        // the interface
        uwall = -nWall & Uwall ;

        //Capillary number fields
        Ca = muwall * uwall / sigmap.value();
    }
    else
    {
        //Info << "Interpolating the PLIC" << nl;
        const volVectorField& interfaceCentre = 
            db.lookupObject<volVectorField>(centresName);

        // Create object for interpolating velocity to isoface centres
        interpolationCellPoint<vector> UInterp(U);

        //Interpolate the cell-centred velocity to the PLIC-centre
        forAll(interfaceCentre, cellI)
        {
            Uplic[cellI] = UInterp.interpolate(interfaceCentre[cellI], cellI);
        }
        //Get patch fields for PLIC velocity
        const auto& pUplic = Uplic.boundaryField()[patchIndex];

        // Calculated the component of the velocity parallel to the wall
        Uwall = pUplic.patchInternalField() - pUplic;

        Uwall -= (nf & Uwall) * nf;

        // Find the direction of the interface parallel to the wall
        nWall = nHat - (nf & nHat) * nf;

        // Normalise nWall
        nWall /= (mag(nWall) + SMALL);

        // Calculate Uwall resolved normal to the interface parallel to
        // the interface
        uwall = -nWall & Uwall;

        //Capillary number fields
        Ca = muwall * uwall / sigmap.value();
    }

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref();   

    // Visualization.
    auto& clangleBoundaryField = contactLineAngle_.boundaryFieldRef();
    auto& clanglePatchField = clangleBoundaryField[patchIndex];
    clanglePatchField *= 0.0;
    contactLineAngle_*=0.0;
        
    // Visualization of the contact angle
    const fvMesh& mesh = nu.mesh(); 
    const auto& faceOwner = mesh.faceOwner();   

    // Maximum contact line capillary number for post-processing
    scalar CaMax= 0;
    // Stabilization technique applied to filter the wisp contribution 
    // If the contact angle changes more than a (hard-coded) user-defined
    // value, changed the contact angle of the particular face with the average value 
    // of the calculated dynamic contact angles
    scalar thetaAvg = 0.0;
    scalar nCLcells = 0.0;

    forAll(thetaf, faceI)
    {   
        if (hasContactLine(faceI))
        {
            if ((alphaPatchField[faceI] > 1e-8) || (alphaPatchField[faceI] <  (1- 1e-8)) )
            {
            
                contactLineAngle_[faceOwner[patch.start() + faceI]] = thetaf[faceI];

                if (mag(Ca[faceI])> CaMax)
                {
                    CaMax = mag(Ca[faceI]);
                }

                //For 1st couple of time stamps, thetaD = theta_initialized
                // if this is not done, it will take thetad=theta0 degrees as the Ca is zero
                // initially which is not correct.
                if ((mag(Ca[faceI])==0) &   (mesh.time().timeIndex() <=1 ))
                {
                    thetaf[faceI] = thetaf[faceI];
                }
                else // apply Kistler's model
                {
                    thetaf[faceI] = min(
                                180 / constant::mathematical::pi * kistler_model(theta0_, Ca[faceI]),
                                scalar(165)
                                ); 
                }
                //Pout << "thetaf in calculation  :" <<  thetaf[faceI] <<nl;
                thetaAvg+=thetaf[faceI];
                nCLcells++;
                //Pout << "thetaf : " << thetaf[faceI] << " Ca: " << mag(Ca[faceI]) << " beta Contri " << beta_ *mag(Ca[faceI]) << " theta m " << thetam << nl;

            }
            clanglePatchField[faceI] = thetaf[faceI];

        }
       // Pout << "thetaf outside calculation  :" <<  thetaf[faceI] <<nl;
    }
    


    // Stabilization of the huge increase in contact angle 
    // due to couple of cells with high capillary number
    if(mesh.time().timeIndex() > 100) 
    {
        scalar globalTheta = returnReduce(thetaAvg, sumOp<scalar>());
        scalar globalnCLCells = returnReduce(nCLcells, sumOp<scalar>());
        if(globalnCLCells!=0)
        {
            scalar gThetaAvg = globalTheta / globalnCLCells;
            //Info << "thetaAvg " << gThetaAvg << nl;
            forAll(thetaf, faceI)
            {
                // 5 is an arbitrary value, if the difference of theta from avg is greater than 20 then rewrite theta to a avg value
                // Avoids the large CA near the equilibirum state
                // if(thetaf[faceI]!=90)
                //     Pout << "Before the cal: " << thetaf[faceI] << " - " << gThetaAvg << " = " << mag(thetaf[faceI] - gThetaAvg) << nl;
                if( (mag(thetaf[faceI] - gThetaAvg) > 5) & (thetaf[faceI] !=90)) 
                {
                    // if(thetaf[faceI]!=90)
                    //    Pout << "went inside " << nl;
                    thetaf[faceI] = gThetaAvg;
                }
                // if(thetaf[faceI]!=90)
                //     Pout << "After the cal: " << thetaf[faceI] << nl;

            }
        }
    }

    //For plotting of contact line Capillary number 
    scalar nonZeroCount = 0;
    // Check if local scalar is non-zero
    if (CaMax != 0)
    {
        nonZeroCount++;
    }
    //Pout << "Individual Max Capillary number: " <<  CaMax << nl;
    if (returnReduce(nonZeroCount, sumOp<scalar>()) !=0)
    {
        Info << "Max Contact Line Capilary number: " << returnReduce(CaMax, sumOp<scalar>())  / returnReduce(nonZeroCount, sumOp<scalar>())<< nl;
    }

    return thetafTmp;
}

void Foam::kistler::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeEntry("theta0", theta0_);
    writeEntry("value", os);
    //fvPatchField<scalar>::writeValueEntry(os);
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
