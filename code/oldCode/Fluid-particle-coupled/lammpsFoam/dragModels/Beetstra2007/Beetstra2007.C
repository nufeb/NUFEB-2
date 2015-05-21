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

#include "Beetstra2007.H"
#include "scalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Beetstra2007, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        Beetstra2007,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Beetstra2007::Beetstra2007
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

    nuf_ = Dnuf_.value();
    rhof_ = Drhof_.value();

    Info<< "Fluid properties dictionary look-up from Drag model: " << endl
        << "Kinematic Viscosity: " << Dnuf_ << "; " << "Density: "
        << Drhof_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Beetstra2007::~Beetstra2007()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//#define DEBUG_JDBEETSTRA

Foam::tmp<Foam::scalarField> Foam::Beetstra2007::Jd
(
    const scalarField& Ur
) const
{

    if
    (
        Ur.size() != alpha_.size()
     || Ur.size() != pd_.size()
    )
        FatalErrorIn("Beetstra2007::Jd() ")
            << "Inconsistent Ur/Alpha/pd."
            << " Ur size: " << Ur.size()
            << " Alpha size: " << alpha_.size()
            << " pd size: " << pd_.size()
            << abort(FatalError) << endl;

    scalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));
/*    scalarField bp = pow(beta, -2.65);
    scalarField Re = max(beta*Ur*pd_/nuf_, scalar(1.0e-3));
    scalarField Cds = 24.0*(1.0 + 0.15*pow(Re, 0.687))/Re;

    forAll(Re, particleI)
    {
        if (Re[particleI] > 1000.0)
        {
            Cds[particleI] = 0.44;
        }
    }

    // Wen and Yu (1966)
    tmp<scalarField> tKWenYu = 0.75*Cds*rhof_*Ur*bp/pd_;
    scalarField& KWenYu = tKWenYu();

*/

    // Determine averaged diameter over all particles
//    scalar rTotD = 0.0;
//    scalar Xi_ = 0.5;
//    X1 = 0.6614;
//    X2 = 0.3386;

    
/*
    int J = 0;
    forAll(beta, speciesJ)
    {
//            rAveD_ = rAveD + Xi_[speciesJ]/pd_[speciesJ];
           rTotD = rTotD + pd_[speciesJ];
           J = J + 1;
    }
    scalar aveD_ = rTotD/J;
  */  

//    scalar aveD_ = 1.0/rAveD;

//    scalarField aveRe = max(beta*Ur*aveD_/nuf_, scalar(1.0e-3)); //in ErgunWenYu not sure if correct
    scalar aveD_ = 0.0022;
    scalarField aveRe = max(Ur*aveD_/nuf_, scalar(1.0e-3)); //TODO porosity beta is missing in formula. Nope Ur is already superficial velocity (Ur = Ua - Ub)

//    scalarField Fdimless = 18.0*nuf_*rhof_/(aveD_*aveD_)*aveRe/aveRe;
//    tmp<scalarField> Fdimless = 18*rhof_*nuf_/sqr(pd_);

    tmp<scalarField> Fdimless = 10.0*alpha_/sqr(beta) + sqr(beta)*(1.0 + 1.5*sqrt(1.0-beta)) + 0.413*aveRe/(24.0*sqr(beta))*((1.0/beta + 3.0*beta*(1.0-beta)+8.4*pow(aveRe,-0.343))/(1.0 + pow(10.0,3.0*alpha_)*pow(aveRe,-0.5-2.0*alpha_)))*18.0*nuf_*rhof_/(aveD_*aveD_);

//    tmp<scalarField> tInterCoeff = Fdimless;
    scalarField& KFdimless = Fdimless();

    forAll(beta, particleJ)
    {
            KFdimless[particleJ] =
//                (beta*pd_[particleJ]/aveD_ + (1.0 - beta)*sqr(pd_[particleJ]/aveD_))*Fdimless*18.0*nuf_*rhof_/aveD_;
                (beta[particleJ]*pd_[particleJ]/aveD_ + (1.0 - beta[particleJ])*sqr(pd_[particleJ]/aveD_) + 0.064*beta[particleJ]*pow(pd_[particleJ]/aveD_,3))*Fdimless()[particleJ];
    }

//    scalarField interCoeff = Fdimless*18*nuf_*rhof_/aveD_;

#ifdef DEBUG_JDBEETSTRA
    Info<< " ==Report====> " << "Beetstra2007::Jd()  " << endl;
    Info<< " alpha: " << alpha_
        << " aveRe: " << aveRe
        << " Ur:" << Ur
        << " pd size: " << pd_.size()
        << " Fdimless: " << Fdimless
        << " Jd: " << tInterCoeff << endl;
#endif

//    return tInterCoeff;
    return Fdimless;
//        return 18*rhof_*nuf_/sqr(pd_);

}

// ************************************************************************* /
