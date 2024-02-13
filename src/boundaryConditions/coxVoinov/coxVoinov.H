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
    Foam::coxVoinov

Description
    A dynamic alphaContactAngle scalar boundary condition
    employs the empirical coxVoinov model for dynamic contact angle 
    ((taken from Dirk Drunding and Anja Lippert dissertation)) i.e.,
    
    theta_d = f_Hoff (Ca + f_Hoff^(-1) (theta_e));
    f_Hoff(x) = arccos(1-2*tanh(5.16(Ca / (1 + 1.31*Ca^0.99))^0.706))

Developed by:
    Muhammad Hassan Asghar
    Mathematical Modeling and Analysis
    TU  Darmstadt


SourceFiles
    coxVoinov.C

\*---------------------------------------------------------------------------*/

#ifndef coxVoinov_H
#define coxVoinov_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class dynamicAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class coxVoinov
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;

        // scaling factor in Cox Voinov 
        // ln(x/L) -> x is the macro length scale and L is the micro length scale
        scalar beta_;


public:

    //- Runtime type information
    TypeName("coxVoinovDCA");


    // Constructors

        //- Construct from patch and internal field
        coxVoinov
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coxVoinov
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coxVoinov
        //  onto a new patch
        coxVoinov
        (
            const coxVoinov&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        coxVoinov
        (
            const coxVoinov&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coxVoinov(*this)
            );
        }

        //- Construct as copy setting internal field reference
        coxVoinov
        (
            const coxVoinov&,
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
                new coxVoinov(*this, iF)
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
