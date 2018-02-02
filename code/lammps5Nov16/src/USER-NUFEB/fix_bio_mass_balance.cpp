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
}

void FixMassBalance::init()
{
  // register fix kinetics with this class
  co2_carbon = pre_co2_carbon = 0;
  glu_carbon = pre_glu_carbon = 0;

  nh3_nitrogen = pre_nh3_nitrogen = 0;
  no2_nitrogen = pre_no2_nitrogen = 0;
  no3_nitrogen = pre_no3_nitrogen = 0;

  total_bmass = pre_total_bmass = 0;

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

  // get biomass concentration (mol/L)
  for (int i = 0; i < nlocal; i++) {
    int pos = kinetics->position(i);
    double rmassCellVol = atom->rmass[i] / vol;
    rmassCellVol /= 24.6;

    if (pos != -1) total_bmass += rmassCellVol;
  }

  if (cflag == 1) c_element_check();
  if (nflag == 1) n_element_check();
  if (mflag == 1) mass_check();

  pre_total_bmass = total_bmass;
  total_bmass = 0;
}

/* ----------------------------------------------------------------------
 mass balance check for Carbon
 ------------------------------------------------------------------------- */

void FixMassBalance::c_element_check() {
//  for (int i = 0; i < kinetics->bgrids; i++) {
//    for (int nu = 1; nu <= nnus; nu++) {
//      if (strcmp(bio->nuName[nu], "co2") == 0) {
//        co2_carbon += nuS[nu][i];
//        //co2_carbon += nuS[nu][i];
//      } else if (strcmp(bio->nuName[nu], "glu") == 0) {
//        glu_carbon += nuS[nu][i];
//      }
//    }
//  }
//
//  pre_co2_carbon = co2_carbon;
//  pre_glu_carbon = glu_carbon;
//
//  co2_carbon = 0;
//  glu_carbon = 0;
}

/* ----------------------------------------------------------------------
 mass balance check for Nitrogen
 ------------------------------------------------------------------------- */

void FixMassBalance::n_element_check() {

  int bgrids = kinetics->bgrids;

  for (int i = 0; i < kinetics->bgrids; i++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuName[nu], "no2") == 0) {
        no2_nitrogen += kinetics->nuS[nu][i];
      } else if (strcmp(bio->nuName[nu], "nh3") == 0) {
        nh3_nitrogen += kinetics->nuS[nu][i];
      }
    }
  }

  double diff_mass = ((total_bmass - pre_total_bmass) * 0.2) / bgrids;
  double diff_no2 = (no2_nitrogen - pre_no2_nitrogen) / bgrids;
  double diff_nh3 = (nh3_nitrogen - pre_nh3_nitrogen) / bgrids;

  double left = fabs(diff_mass) + fabs(diff_no2);
  double right = fabs(diff_nh3);

  printf("(N) Diff = %e, Biomass = %e, NO2 = %e, NH3  = %e \n",
      left-right, diff_mass, diff_no2, diff_nh3);

  pre_nh3_nitrogen = nh3_nitrogen;
  pre_no2_nitrogen = no2_nitrogen;

  nh3_nitrogen = 0;
  no2_nitrogen = 0;
}

void FixMassBalance::mass_check() {

}
