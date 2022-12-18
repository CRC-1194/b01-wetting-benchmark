/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "mysetFlow.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mysetFlow, 0);
    addToRunTimeSelectionTable(functionObject, mysetFlow, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::mysetFlow::modeType
>
Foam::functionObjects::mysetFlow::modeTypeNames
({
    { functionObjects::mysetFlow::modeType::MYVORTEX2D, "myvortex2D" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mysetFlow::setPhi(const volVectorField& U)
{
    surfaceScalarField* phiptr =
        mesh_.getObjectPtr<surfaceScalarField>(phiName_);

    if (!phiptr)
    {
        return;
    }

    if (rhoName_ != "none")
    {
        const volScalarField* rhoptr =
            mesh_.findObject<volScalarField>(rhoName_);

        if (rhoptr)
        {
            const volScalarField& rho = *rhoptr;
            *phiptr = fvc::flux(rho*U);
        }
        else
        {
            FatalErrorInFunction
                << "Unable to find rho field'" << rhoName_
                << "' in the mesh database.  Available fields are:"
                << mesh_.names<volScalarField>()
                << exit(FatalError);
        }
    }
    else
    {
        *phiptr = fvc::flux(U);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mysetFlow::mysetFlow
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mode_(modeType::MYVORTEX2D),
    UName_("U"),
    rhoName_("none"),
    phiName_("phi"),
    reverseTime_(VGREAT),
    scalePtr_(nullptr),
    origin_(Zero),
    R_(tensor::I),
    omegaPtr_(nullptr),
    velocityPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mysetFlow::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Info<< name() << ":" << endl;

        modeTypeNames.readEntry("mode", dict, mode_);

        Info<< "    operating mode: " << modeTypeNames[mode_] << endl;

        if (dict.readIfPresent("U", UName_))
        {
            Info<< "    U field name: " << UName_ << endl;
        }

        if (dict.readIfPresent("rho", rhoName_))
        {
            Info<< "    rho field name: " << rhoName_ << endl;
        }

        if (dict.readIfPresent("phi", phiName_))
        {
            Info<< "    phi field name: " << phiName_ << endl;
        }

        if (dict.readIfPresent("reverseTime", reverseTime_))
        {
            Info<< "    reverse flow direction at time: " << reverseTime_
                << endl;
            reverseTime_ = mesh_.time().userTimeToTime(reverseTime_);
        }

        // Scaling is applied across all modes
        scalePtr_ = Function1<scalar>::New("scale", dict);

        switch (mode_)
        {
            
            case modeType::MYVORTEX2D:
            {
                dict.readEntry("origin", origin_);
                const vector refDir(dict.get<vector>("refDir").normalise());
                const vector axis(dict.get<vector>("axis").normalise());

                R_ = tensor(refDir, axis, refDir^axis);
                break;
            }
        }

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::mysetFlow::execute()
{
    volVectorField* Uptr =
        mesh_.getObjectPtr<volVectorField>(UName_);

    surfaceScalarField* phiptr =
        mesh_.getObjectPtr<surfaceScalarField>(phiName_);

    Log << nl << name() << ":" << nl;

    if (!Uptr || !phiptr)
    {
        Info<< "    Either field " << UName_ << " or " << phiName_
            << " not found in the mesh database" << nl;

        return true;
    }

    const scalar t = mesh_.time().timeOutputValue();

    Log << "    setting " << UName_ << " and " << phiName_ << nl;

    volVectorField& U = *Uptr;
    surfaceScalarField& phi = *phiptr;

    switch (mode_)
    {
        case modeType::MYVORTEX2D:
        {
            const scalar pi = Foam::constant::mathematical::pi;

            const volVectorField& C = mesh_.C();

            const volVectorField d
            (
                typeName + ":d",
                C - dimensionedVector("origin", dimLength, origin_)
            );
            const scalarField x(d.component(vector::X));
            const scalarField z(d.component(vector::Z));

            //Vortex-in-a-box velocity field 
            // prescribed for a 2D simulation
            vectorField& Uc = U.primitiveFieldRef();
            Uc.replace(vector::X, -sin(pi*x)*cos(pi*z));
            Uc.replace(vector::Y, scalar(0));
            Uc.replace(vector::Z, cos(pi*x)*sin(pi*z));

            U = U & R_;

            // Calculating incompressible flux based on stream function
            // Note: R_ rotation not implemented in phi calculation
            const scalarField xp(mesh_.points().component(0) - origin_[0]);
            const scalarField yp(mesh_.points().component(1) - origin_[1]);
            const scalarField zp(mesh_.points().component(2) - origin_[2]);
            const scalarField psi((1.0/pi)*sqr(sin(pi*xp))*sqr(sin(pi*zp)));

            scalarField& phic = phi.primitiveFieldRef();
            forAll(phic, fi)
            {
                phic[fi] = 0;
                const face& f = mesh_.faces()[fi];
                const label nPoints = f.size();

                forAll(f, fpi)
                {
                    const label p1 = f[fpi];
                    const label p2 = f[(fpi + 1) % nPoints];
                    phic[fi] += 0.5*(psi[p1] + psi[p2])*(yp[p2] - yp[p1]);
                }
            }

            surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();
            forAll(phibf, patchi)
            {
                scalarField& phif = phibf[patchi];
                const label start = mesh_.boundaryMesh()[patchi].start();

                forAll(phif, fi)
                {
                    phif[fi] = 0;
                    const face& f = mesh_.faces()[start + fi];
                    const label nPoints = f.size();

                    forAll(f, fpi)
                    {
                        const label p1 = f[fpi];
                        const label p2 = f[(fpi + 1) % nPoints];
                        phif[fi] += 0.5*(psi[p1] + psi[p2])*(yp[p2] - yp[p1]);
                    }
                }
            }

            break;
        }
    }

    if (t > reverseTime_)
    {
        Log << "    flow direction: reverse" << nl;
        U.negate();
        phi.negate();
    }

    // Apply scaling
    const scalar s = scalePtr_->value(t);
    U *= s;
    phi *= s;

    U.correctBoundaryConditions();

    const scalarField sumPhi(fvc::surfaceIntegrate(phi));
    Log << "    Continuity error: max(mag(sum(phi))) = "
        << gMax(mag(sumPhi)) << nl << endl;

    return true;
}


bool Foam::functionObjects::mysetFlow::write()
{
    const auto* Uptr = mesh_.findObject<volVectorField>(UName_);
    if (Uptr)
    {
        Uptr->write();
    }

    const auto* phiptr = mesh_.findObject<surfaceScalarField>(phiName_);
    if (phiptr)
    {
        phiptr->write();
    }

    return true;
}


// ************************************************************************* //
