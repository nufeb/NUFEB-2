/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_bio_mass_balance.h"

#include <cstdio>
#include <cstring>

#include "atom.h"
#include "bio.h"
#include "error.h"
#include "fix_bio_kinetics.h"
#include "force.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include <stdio.h>
#include <math.h>


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMassBalance::FixMassBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg > 7) error->all(FLERR,"Illegal fix mbalance command");

  nevery = force->inumeric(FLERR,arg[3]);

  int nkeywords = narg - 4;
  nflag = cflag = mflag = 0;

  for (int i = 0; i < nkeywords; i++) {
    if (strcmp(arg[4+i], "n") == 0) nflag = 1;
    else if (strcmp(arg[4+i], "c") == 0) cflag = 1;
    else if (strcmp(arg[4+i], "mass") == 0) mflag = 1;
    else error->all(FLERR,"Illegal fix mbalance keywords");
  }

}

FixMassBalance::~FixMassBalance()
{
  delete bmass;
  delete pre_bmass;
}

void FixMassBalance::init()
{
  // register fix kinetics with this class
  co2_carbon = pre_co2_carbon = 0;
  nh3mass = pre_nh3mass = 0;
  glu_carbon = pre_glu_carbon = 0;
  kinetics = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for mass balance check");

  bio = kinetics->bio;
  vol = kinetics->stepx * kinetics->stepy * kinetics->stepz;

  bmass = new double[kinetics->bgrids]();
  pre_bmass = new double[kinetics->bgrids]();
}

/* ---------------------------------------------------------------------- */

int FixMassBalance::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


void FixMassBalance::end_of_step() {
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;

  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;

  gYield = kinetics->gYield;
  nuS = kinetics->nuS;
  nnus = bio->nnus;

  // get biomass concentration at each grid
  for (int i = 0; i < nlocal; i++) {
      int pos = kinetics->position(i);
      double rmassCellVol = atom->rmass[i] / vol;
      rmassCellVol /= 24.6;

      if (pos != -1) bmass[pos] += rmassCellVol;
  }

  if (cflag == 1) c_element_check();
  if (nflag == 1) n_element_check();
  if (mflag == 1) mass_check();

  for (int i = 0; i < kinetics->bgrids; i++) {
    pre_bmass[i] = bmass[i];
    bmass[i] = 0;
  }
}

void FixMassBalance::c_element_check() {

  double tbmass = 0;  //total biomass
  double pre_tbmass = 0;


  for (int i = 0; i < kinetics->bgrids; i++) {

    tbmass += bmass[i];
    pre_tbmass += pre_bmass[i];

    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuName[nu], "co2") == 0) {
        co2_carbon += nuS[nu][i];
        //co2_carbon += nuS[nu][i];
      } else if (strcmp(bio->nuName[nu], "glu") == 0) {
        glu_carbon += nuS[nu][i];
      }
    }
 //   double* invYield = new double[atom->ntypes+1];

 //   for (int j = 1; j < atom->ntypes+1; j++) {
//      if (gYield[j][i] != 0) invYield[j] = 1 / gYield[j][i];
//      else invYield = 0;
//
//      for (int nu = 1; nu <= nnus; nu++) {
//        if (strcmp(bio->nuName[nu], "co2") == 0) {
//          // yield of co2
//          double metCoeff = catCoeff[j][nu] * invYield[j] + anabCoeff[j][nu];
//          // weight of co2
//          co2_carbon += nuS[nu][i];
//          //co2_carbon += nuS[nu][i];
//        } else if (strcmp(bio->nuName[nu], "glu") == 0) {
//          double metCoeff = catCoeff[j][nu] * invYield[j] + anabCoeff[j][nu];
//          // weight of glu
//          //glu_carbon += nuS[nu][i] * 6;
//          glu_carbon += nuS[nu][i] * 6;
//        }
//      }
//    }


//    delete invYield;
  }

  double left = fabs((tbmass - pre_tbmass) ) + fabs((co2_carbon - pre_co2_carbon));
  double right = fabs((glu_carbon - pre_glu_carbon) * 6);

  printf("mass element diff is = %e, mass = %e, co2 = %e, glu  = %e \n",
      left-right, fabs((tbmass - pre_tbmass)), fabs((co2_carbon - pre_co2_carbon)), fabs((glu_carbon - pre_glu_carbon)));

  pre_co2_carbon = co2_carbon;
  pre_glu_carbon = glu_carbon;

  co2_carbon = 0;
  glu_carbon = 0;
}

void FixMassBalance::n_element_check() {

}

void FixMassBalance::mass_check() {

  for (int i = 0; i < kinetics->bgrids; i++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuName[nu], "co2") == 0) {
        // weight of co2
        co2_carbon += nuS[nu][i] * vol * 1000;
      } else if (strcmp(bio->nuName[nu], "glu") == 0) {
        glu_carbon += nuS[nu][i] * vol * 1000;
      } else if (strcmp(bio->nuName[nu], "nh3") == 0) {
        nh3mass += nuS[nu][i] * vol * 1000;
      } else if (strcmp(bio->nuName[nu], "o2") == 0) {
        o2mass += nuS[nu][i] * vol * 1000;
      }
    }
  }

  // unit for mol to kg
  glu_carbon = glu_carbon * 180.156 / 1000;
  nh3mass = nh3mass * 17.031 / 1000;
  o2mass = o2mass * 15.999 / 1000;
  co2_carbon = co2_carbon * 44.01 / 1000;


  double left = (bmass - pre_bmass) * 24.6 / 1000;
  double right = (o2mass - pre_o2mass) + (glu_carbon - pre_glu_carbon) + (nh3mass - pre_nh3mass) + (co2_carbon - pre_co2_carbon);

  printf("mass element diff is = %e, right = %e, left = %e \n", left-right, right, left);

  pre_co2_carbon = co2_carbon;
  pre_glu_carbon = glu_carbon;
  pre_nh3mass = nh3mass;
  pre_o2mass = o2mass;

  co2_carbon = 0;
  glu_carbon = 0;
  nh3mass = 0;
  o2mass = 0;
}
