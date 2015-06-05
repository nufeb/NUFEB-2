/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Stokes.H"
#include "scalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Stokes, 0);

    addToRunTimeSelectionTable
    (
        dragModel, 
        Stokes,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Stokes::Stokes
(
    const dictionary& cloudDict,
    const IOdictionary& transDict,
    const scalarField& alpha,
    const scalarField& pd
)
:
    dragModel
    (
          cloudDict,
          transDict,
          alpha,
          pd
    )
{
    dimensionedScalar Dnuf_(transDict_.lookup("nub"));
    dimensionedScalar Drhof_(transDict_.lookup("rhob"));
//    dimensionedScalar Dda_(transDict_.lookup("da"));

    nuf_ = Dnuf_.value();
    rhof_ = Drhof_.value();
//    da_ = Dda_.value();

    Info<< "Fluid properties dictionary look-up from Drag model: " << endl
        << "Kinematic Viscosity: " << Dnuf_ << "; " << "Density: "
        << Drhof_ << endl;//"; " << "Particle diameter: " << Dda_ << 
//        endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Stokes::~Stokes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//#define DEBUG_JDSTOKES

Foam::tmp<Foam::scalarField> Foam::Stokes::Jd
(
    const scalarField& Ur
) const
{

    if
    (
        Ur.size() != alpha_.size()
     || Ur.size() != pd_.size()
    )
    {
        FatalErrorIn("Stokes::Jd() ")
            << "Inconsistent Ur/Alpha/pd."
            << " Ur size: " << Ur.size()
            << " Alpha size: " << alpha_.size()
            << " pd size: " << pd_.size()
            << abort(FatalError) << endl;
    }

    
    scalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));

    scalarField Re = max(Ur*pd_/nuf_, scalar(1.0e-3));

    scalarField Cds = 24.0/Re + 2.6*(Re/5.0)/(1.0+pow(Re/5.0,1.52)) + 0.411*pow(Re/263000,-7.94)/(1+pow(Re/263000,-8.0)) + (pow(Re,0.8)/461000);

#ifdef DEBUG_JDSTOKES
    Info<< " ==Report====> " << "Stokes::Jd()  " << endl;
    Info<< "Lagrangian rho: " << rhof_
        << "Lagrangian nu: " << nuf_
        << "Lagrangian pd: " << pd_
        << "Lagrangian Re: " << Re
        << "Lagrangian Ur:" << Ur
        << "Lagrangian Cds: " << Cds << endl;
//        << "Lagrangian Jd: " <<  0.5*Cds*rhof_*sqr(Ur)*(constant::mathematical::pi*sqr(pd_/2)) << endl;
#endif

        return 18*rhof_*nuf_/sqr(pd_);
//          return 0.75*rhof_/pd_*Cds*Ur;
}

// ************************************************************************* //
