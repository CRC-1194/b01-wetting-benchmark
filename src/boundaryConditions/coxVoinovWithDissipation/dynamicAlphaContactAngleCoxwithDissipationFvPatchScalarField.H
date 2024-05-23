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

Class
    Foam::dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField

Description
    A dynamic alphaContactAngle scalar boundary condition
    employs the hydrodynamic CoxwithDissipation model for dynamic contact angle 
    with uncompensated Young-stress term for dissipation
    

    ((taken from Dirk Drunding and Anja Lippert dissertation)) i.e.,
    
    // thetad = (beta*  Ca + theta0^3)^(1/3) for \theta<135° (Dirk Grunding dissertation)
    // beta = ln(x/L) -> x is the macro length scale and L is the micro length scale

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt


SourceFiles
    dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField_H
#define dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "volFields.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class CoxwithDissipation Declaration
\*---------------------------------------------------------------------------*/

class dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;
        scalar thetaA_;
        scalar thetaR_;
        scalar diss_;

        // scaling factor in Cox Voinov 
        // ln(x/L) -> x is the macro length scale and L is the micro length scale
        scalar beta_;

        //- Visualization 
        mutable volScalarField contactLineAngle_;

    // Private member functions
       
        //- Determine whether a face contains a contact line
        bool hasContactLine
        (
            label faceI
        ) const;


public:

    //- Runtime type information
    TypeName("dynamicAlphaContactAngleCoxwithDissipation");


    // Constructors

        //- Construct from patch and internal field
        dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  CoxwithDissipation
        //  onto a new patch
        dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField
        (
            const dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField
        (
            const dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField
        (
            const dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicAlphaContactAngleCoxwithDissipationFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate and return dynamic contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
