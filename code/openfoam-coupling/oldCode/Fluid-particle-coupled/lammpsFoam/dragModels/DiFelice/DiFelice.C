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

#include "DiFelice.H"
#include "scalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DiFelice, 0);

    addToRunTimeSelectionTable
    (
        dragModel, 
        DiFelice,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DiFelice::DiFelice
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

Foam::DiFelice::~DiFelice()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//#define DEBUG_DIFELICE

Foam::tmp<Foam::scalarField> Foam::DiFelice::Jd
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
        FatalErrorIn("DiFelice::Jd() ")
            << "Inconsistent Ur/Alpha/pd."
            << " Ur size: " << Ur.size()
            << " Alpha size: " << alpha_.size()
            << " pd size: " << pd_.size()
            << abort(FatalError) << endl;
    }

    
    scalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));

    scalarField Re = max(Ur*pd_/nuf_, scalar(1.0e-3));

    scalarField Cd0 = sqr((0.63 + 4.8/pow(Re,0.5)));
    scalarField voidF = pow(beta,(-(3.7-0.65*exp(-0.5*sqr(1.5-log10(Re))))));



#ifdef DEBUG_DIFELICE
    Info<< " ==Report====> " << "Stokes::Jd()  " << endl;
    Info<< "Lagrangian rho: " << rhof_
        << "Lagrangian nu: " << nuf_
        << "Lagrangian pd: " << pd_
        << "Lagrangian Re: " << Re
        << "Lagrangian Ur:" << Ur
        << "Lagrangian Cds: " << Cd0 // << endl;
        << "voidF: " << voidF
        << "Lagrangian Jd: " << 1.5*rhof_/pow(pd_,3)*Ur*Cd0*voidF << endl;
#endif

        return 1.5*rhof_/pd_*Ur*Cd0*voidF;

}

// ************************************************************************* //
