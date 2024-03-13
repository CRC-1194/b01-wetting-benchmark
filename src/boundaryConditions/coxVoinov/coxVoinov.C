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

#include "interpolationCellPoint.H"

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

bool Foam::coxVoinov::
hasContactLine
(
    label faceI
) const
{
        // Look up PLIC normals and positions. 
    const auto& db = this->db(); 

    const auto normalsName = IOobject::groupName
    (
        "interfaceNormal", 
        this->internalField().group()
    );
    const auto centresName = IOobject::groupName
    (
        "interfaceCentre", 
        this->internalField().group()
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

Foam::coxVoinov::
coxVoinov
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    theta0_(0.0),
    beta_(0.0),
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
    beta_(gcpsf.beta_),
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
    beta_(readScalar(dict.lookup("beta"))),
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


Foam::coxVoinov::
coxVoinov
(
    const coxVoinov& gcpsf
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf),
    theta0_(gcpsf.theta0_),
    beta_(gcpsf.beta_),
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


Foam::coxVoinov::
coxVoinov
(
    const coxVoinov& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    beta_(gcpsf.beta_),
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

Foam::tmp<Foam::scalarField>
Foam::coxVoinov::theta
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

    const fvPatchScalarField & nu1p = nu1.boundaryField()[patchIndex];

    word fieldName = "alpha." + phase1Name;
    const volScalarField& alpha1_ =
        nu1.mesh().lookupObject<volScalarField>(fieldName);

    const volScalarField::Boundary& abf = alpha1_.boundaryField();   
    // Lookup the desired alpha values on the patch
    const fvPatchField<scalar>& alphaPatchField = abf[patchIndex];
    
    //Initialization without velocity which is updated later in the code
    //Capillary number fields
    scalarField Ca(nu1p * rho1.value() / sigmap.value());

    // 1. Get the PLIC centre and normals
    // Look up PLIC normals and positions. 
    const auto& db = this->db(); 

    const auto normalsName = IOobject::groupName
    (
        "interfaceNormal", 
        this->internalField().group()
    );
    const auto centresName = IOobject::groupName
    (
        "interfaceCentre", 
        this->internalField().group()
    );

    bool hasNormals = db.found(normalsName);
    bool hasCentres = db.found(centresName);

    //If there are no PLIC elements i.e, no wall velocity making Capillary number equals zero
    // Else calculate Ca using interpolated cell centred velocity to PLIC centre (using inverse distance weightings)
    if (!hasNormals || !hasCentres )
    {
        // This BC is updated before interface reconstruction.
        // Do nothing if PLIC fields are not available in the registry. 

        // Calculated the component of the velocity parallel to the wall
        vectorField Uwall(Up.patchInternalField() - Up);
        //Info << "Patch Internal Field: \n " << Uwall<< nl;

        //Info << "\n Patch internal field: \n" << Up.patchInternalField();
        Uwall -= (nf & Uwall) * nf;

        // Find the direction of the interface parallel to the wall
        vectorField nWall(nHat - (nf & nHat) * nf);

        // Normalise nWall
        nWall /= (mag(nWall) + SMALL);

        // Calculate Uwall resolved normal to the interface parallel to
        // the interface
        scalarField uwall(-nWall & Uwall);

        //Capillary number fields
        Ca = nu1p * rho1.value() * uwall / sigmap.value();
    }
    else
    {
        const volVectorField& interfaceCentre = 
            db.lookupObject<volVectorField>(centresName);

        // Cell-centred velocity field
        Foam::volVectorField& U  = const_cast<Foam::volVectorField&> (this -> db().objectRegistry::lookupObject < volVectorField > ("U"));

        // velocity field at the plic centre
        // will be updated using inverse distance interpolation scheme on U
        //This field is then used to calculate the Capillary number as in isoAdvection.C
        // plicU is used for advection
        Foam::volVectorField Uplic (U*scalar(0));

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
        vectorField Uwall(pUplic.patchInternalField() - pUplic);
        //Info << "Patch Internal Field: \n " << Uwall<< nl;

        //Info << "\n Patch internal field: \n" << Up.patchInternalField();
        Uwall -= (nf & Uwall) * nf;

        // Find the direction of the interface parallel to the wall
        vectorField nWall(nHat - (nf & nHat) * nf);

        // Normalise nWall
        nWall /= (mag(nWall) + SMALL);

        // Calculate Uwall resolved normal to the interface parallel to
        // the interface
        scalarField uwall(-nWall & Uwall);

        //Capillary number fields
        Ca = nu1p * rho1.value() * uwall / sigmap.value();
    }

    // Compute the contact angles at the wall.
    tmp<scalarField> thetafTmp = Foam::radToDeg(Foam::acos(nHat & nf));
    scalarField& thetaf = thetafTmp.ref();   

    // Visualization.
    auto& clangleBoundaryField = contactLineAngle_.boundaryFieldRef();
    auto& clanglePatchField = clangleBoundaryField[patchIndex];
    clanglePatchField *= 0.0;
        
    // Visualization of the contact angle
    const fvMesh& mesh = nu1.mesh(); 
    const auto& faceOwner = mesh.faceOwner();     

    // For all boundary faces
    forAll(thetaf, faceI)
    {    
        if (hasContactLine(faceI) )
        {
            // Visualize current PLIC contact angle as a cell-centered value.
            contactLineAngle_[faceOwner[patch.start() + faceI]] = thetaf[faceI];

            thetaf[faceI] = min(
            180 / constant::mathematical::pi * (pow(beta_  * mag(Ca[faceI]) +
                pow(theta0_ * constant::mathematical::pi / 180, 3),
                0.3333333)),
            scalar(135)
        );
            clanglePatchField[faceI] = thetaf[faceI];
        }
    }

    return thetafTmp;
    // return min(
    //         180 / constant::mathematical::pi * (pow(beta_ * pos(Ca) * Ca +
    //             pow(theta0_ * constant::mathematical::pi / 180, 3),
    //             0.3333333)),
    //         scalar(135)
    //     );

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
