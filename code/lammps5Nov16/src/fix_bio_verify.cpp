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

#include "fix_bio_verify.h"

#include <cstdio>
#include <cstring>

#include "atom.h"
#include "bio.h"
#include "error.h"
#include "fix_bio_kinetics.h"
#include "compute_bio_height.h"
#include "fix_bio_kinetics_diffusion.h"
#include "fix_bio_kinetics_monod.h"
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

FixVerify::FixVerify(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
 // if (narg > 7) error->all(FLERR,"Illegal fix verify command");

  nevery = force->inumeric(FLERR,arg[3]);

  int nkeywords = narg - 4;
  bm1flag = bm2flag = bm3flag = mflag = 0;
  bm1cflag = 0;

  for (int i = 0; i < nkeywords; i++) {
    if (strcmp(arg[4+i], "bm1") == 0) {
      bm1flag = 1;
      bm1cflag = force->inumeric(FLERR,arg[5+i]);
    }
    else if (strcmp(arg[4+i], "bm2") == 0) bm2flag = 1;
    else if (strcmp(arg[4+i], "mb") == 0) mflag = 1;
    else if (strcmp(arg[4+i], "bm3") == 0) bm3flag = 1;
  }

}

FixVerify::~FixVerify()
{
}

void FixVerify::init()
{
  nh3_nitrogen = pre_nh3_nitrogen = 0;
  no2_nitrogen = pre_no2_nitrogen = 0;
  no3_nitrogen = pre_no3_nitrogen = 0;

  total_bmass = pre_total_bmass = 0;

  kinetics = NULL;

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required");

  bio = kinetics->bio;
  vol = kinetics->stepx * kinetics->stepy * kinetics->stepz;
  kinetics->diffusion->bulkflag = 0;
}

/* ---------------------------------------------------------------------- */

int FixVerify::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


void FixVerify::end_of_step() {
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;

  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;

  gYield = kinetics->gYield;
  nuS = kinetics->nuS;
  nnus = bio->nnus;

  if (mflag == 1) nitrogen_mass_balance();
  if (bm1flag == 1) benchmark_one();
}

/* ----------------------------------------------------------------------
 mass balance check for Nitrogen
 ------------------------------------------------------------------------- */

void FixVerify::nitrogen_mass_balance() {
  // get biomass concentration (mol/L)
  for (int i = 0; i < nlocal; i++) {
    int pos = kinetics->position(i);
    double rmassCellVol = atom->rmass[i] / vol;
    rmassCellVol /= 24.6;

    if (pos != -1) total_bmass += rmassCellVol;
  }

  // get nitrogen concentration
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

  pre_total_bmass = total_bmass;
  total_bmass = 0;
}

/* ----------------------------------------------------------------------
 biofilm benchmark problem one
 ------------------------------------------------------------------------- */

void FixVerify::benchmark_one() {
  // register compute with this class
  class ComputeNufebHeight *cheight;
  int ncompute = modify->ncompute;
  int nfix = modify->nfix;

  for (int j = 0; j < ncompute; j++) {
    if (strcmp(modify->compute[j]->style, "ave_height") == 0) {
      cheight = static_cast<ComputeNufebHeight *>(lmp->modify->compute[j]);
      break;
    }
  }

  // case 3, average biofilm thickness: Lf = 20 μm
  if (bm1cflag == 3) {
    // solve mass balance in bulk liquid
    double ave_height = cheight->compute_scalar();

    if (ave_height > 0.2e-4) {
      kinetics->diffusion->bulkflag = 1;
    }

    // cease growth and division
    if (ave_height > 0.145e-4) {
      kinetics->monod->external_gflag = 0;
    }
  }
}

