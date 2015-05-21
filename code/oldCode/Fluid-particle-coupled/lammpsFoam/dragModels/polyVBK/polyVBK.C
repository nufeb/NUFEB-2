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

#include "polyVBK.H"
#include "scalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyVBK, 0);

    addToRunTimeSelectionTable
    (
        polyDragModel,
        polyVBK,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyVBK::polyVBK
(
    const dictionary& cloudDict,
    const IOdictionary& transDict,
    const scalarField& alpha,
    const scalarField& pd,
    const scalarField& aveD
)
:
    polyDragModel
    (
        cloudDict,
        transDict,
        alpha,
        pd,
        aveD
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

Foam::polyVBK::~polyVBK()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//#define DEBUG_JDBEETSTRA

Foam::tmp<Foam::scalarField> Foam::polyVBK::Jd
(
    const scalarField& Ur, const scalarField& sauterAveD
) const
{

    if
    (
        Ur.size() != alpha_.size()
     || Ur.size() != pd_.size()
    )
        FatalErrorIn("polyVBK::Jd() ")
            << "Inconsistent Ur/Alpha/pd."
            << " Ur size: " << Ur.size()
            << " Alpha size: " << alpha_.size()
            << " pd size: " << pd_.size()
            << abort(FatalError) << endl;

    scalarField beta = max(scalar(1) - alpha_, scalar(1.0e-6));
//   Info<< "Debug 71c" << endl;

//    scalarField aveD_ = sauterAveD;
    scalarField aveRe = max(Ur*aveD_/nuf_, scalar(1.0e-3)); //TODO porosity beta is missing in formula. Nope Ur is already superficial velocity (Ur = Ua - Ub)
//   Info<< "Debug 71d" << endl;



   Info<< "Debug 71X" << endl;
//   scalarField M = aveD_;
/*    forAll(aveD_, cI){
      if(aveD_[cI] < scalar(1.0e-6)){
        M[cI] = 0.0;
      }else{
        M[cI] = 18.0*nuf_*rhof_/sqr(aveD_[cI]);//aveD_);
      }
    }   */
   scalarField aveDmin = max(aveD_, scalar(1.0e-6));
   Info<< "Debug 71XY" << endl;
   scalarField yi = pd_/aveDmin;
   Info<< "Debug 71XYY" << endl;
   scalarField nonDim = 18.0*rhof_*nuf_/sqr(aveDmin);

   Info<< "Debug 71XX" << endl;

   scalarField monoF = 10.0*alpha_/sqr(beta) + sqr(beta)*(1.0 + 1.5*sqrt(1.0-beta)) + 0.413*aveRe/(24.0*sqr(beta))*((1.0/beta + 3.0*beta*(1.0-beta)+8.4*pow(aveRe,-0.343))/(1.0 + pow(10.0,3.0*alpha_)*pow(aveRe,-0.5-2.0*alpha_)));
   scalarField polyFactor = beta*yi + (1.0 - beta)*sqr(yi) + 0.064*beta*pow(yi,3);

   Info<< "Debug 71XXX" << endl;
//   tmp<scalarField> polyFdimless = polyFactor*monoF*nonDim;
   tmp<scalarField> polyFdimless = monoF*nonDim;

/*   tmp<scalarField> Fdimless = 10.0*alpha_/sqr(beta) + sqr(beta)*(1.0 + 1.5*sqrt(1.0-beta)) + 0.413*aveRe/(24.0*sqr(beta))*((1.0/beta + 3.0*beta*(1.0-beta)+8.4*pow(aveRe,-0.343))/(1.0 + pow(10.0,3.0*alpha_)*pow(aveRe,-0.5-2.0*alpha_)))*M;

   Info<< "Debug 71e" << endl;
    scalarField& KFdimless = Fdimless();

    forAll(beta, particleJ)
    {
        if(aveDmin[particleJ] >= scalar(1.0e-6)){
            KFdimless[particleJ] =
                (beta[particleJ]*pd_[particleJ]/aveDmin[particleJ] + (1.0 - beta[particleJ])*sqr(pd_[particleJ]/aveDmin[particleJ]) + 0.064*beta[particleJ]*pow(pd_[particleJ]/aveDmin[particleJ],3))*Fdimless()[particleJ];
        }
    }
*/
   Info<< "Debug 71f" << endl;


#ifdef DEBUG_JDBEETSTRA
    Info<< " ==Report====> " << "polyVBK::Jd()  " << endl;
    Info<< " alpha: " << alpha_
        << " aveRe: " << aveRe
        << " Ur:" << Ur
        << " pd size: " << pd_.size()
	<< " aveDmin: " << aveDmin
        << " Fdimless: " << Fdimless << endl;
//        << " Jd: " << tInterCoeff << endl;
#endif

    return polyFdimless;

}

// ************************************************************************* /
