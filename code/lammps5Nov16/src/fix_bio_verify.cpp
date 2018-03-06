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
#include "atom_vec_bio.h"
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
#include "comm.h"


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
  demflag = 0;

  for (int i = 0; i < nkeywords; i++) {
    if (strcmp(arg[4+i], "bm1") == 0) {
      bm1flag = 1;
      bm1cflag = force->inumeric(FLERR,arg[5+i]);
    }
    else if (strcmp(arg[4+i], "bm2") == 0) bm2flag = 1;
    else if (strcmp(arg[4+i], "mb") == 0) mflag = 1;
    else if (strcmp(arg[4+i], "bm3") == 0) bm3flag = 1;
    else if (strcmp(arg[4+i], "demflag") == 0) {
      if (strcmp(arg[5+i],"no") == 0) demflag = 0;
      else if (strcmp(arg[5+i],"yes") == 0) demflag = 1;
      else error->all(FLERR,"Illegal demflag parameter (yes or no)");
    }
  }
}

FixVerify::~FixVerify()
{
}

void FixVerify::init()
{
  // get overall concentration
  global_no2 = 0;
  global_nh3 = 0;
  global_pre_no2 = 0;
  global_pre_nh3 = 0;
  global_smass = 0;
  global_pre_smass = 0;

  kinetics = NULL;
  diffusion = NULL;

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/diffusion") == 0) {
      diffusion = static_cast<FixKineticsDiffusion *>(lmp->modify->fix[j]);
    }
  }

  int ncompute = modify->ncompute;

  for (int j = 0; j < ncompute; j++) {
    if (strcmp(modify->compute[j]->style, "ave_height") == 0) {
      cheight = static_cast<ComputeNufebHeight *>(lmp->modify->compute[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required");

  bio = kinetics->bio;
  vol = kinetics->stepx * kinetics->stepy * kinetics->stepz;
 // kinetics->diffusion->bulkflag = 0;
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
  if (demflag) return;

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
  double smass, pre_smass;
  double nh3_nitrogen, pre_nh3_nitrogen;
  double no2_nitrogen, pre_no2_nitrogen;
  double no3_nitrogen, pre_no3_nitrogen;
  int ngrids;

  nh3_nitrogen = pre_nh3_nitrogen = 0;
  no2_nitrogen = pre_no2_nitrogen = 0;
  no3_nitrogen = pre_no3_nitrogen = 0;

  smass = 0;

  // get biomass concentration (mol/L)
  for (int i = 0; i < nlocal; i++) {
    int pos = kinetics->position(i);
    double rmassCellVol = atom->rmass[i] / vol;
    rmassCellVol /= 24.6;

    if (pos != -1) smass += rmassCellVol;
  }

  // get overall biamass concentration
  MPI_Allreduce(&smass,&global_smass,1,MPI_DOUBLE,MPI_SUM,world);

  // get nitrogen concentration
  int bgrids = kinetics->bgrids;

  for (int i = 0; i < bgrids; i++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuName[nu], "no2") == 0) {
        no2_nitrogen += kinetics->nuS[nu][i];
      } else if (strcmp(bio->nuName[nu], "nh3") == 0) {
        nh3_nitrogen += kinetics->nuS[nu][i];
      }
    }
  }

  no2_nitrogen /= bgrids;
  nh3_nitrogen /= bgrids;

  MPI_Allreduce(&no2_nitrogen,&global_no2,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&nh3_nitrogen,&global_nh3,1,MPI_DOUBLE,MPI_SUM,world);

  double diff_mass = ((global_smass - global_pre_smass) * 0.2) / (kinetics->nx * kinetics->ny * kinetics->nz);
  double diff_no2 = (global_no2 - global_pre_no2) / comm->nprocs;
  double diff_nh3 = (global_nh3 - global_pre_nh3) / comm->nprocs;

  double left = fabs(diff_mass) + fabs(diff_no2);
  double right = fabs(diff_nh3);

  if (comm->me == 0) printf("(N) Diff = %e, Biomass = %e, NO2 = %e, NH3  = %e \n",
      left-right, diff_mass, diff_no2, diff_nh3);

  global_pre_nh3 = global_nh3;
  global_pre_no2 = global_no2;
  global_pre_smass = global_smass;

  global_nh3 = 0;
  global_no2 = 0;
  global_smass = 0;
}

/* ----------------------------------------------------------------------
 biofilm benchmark problem one
 ------------------------------------------------------------------------- */

void FixVerify::benchmark_one() {
  double global_tmass;
  double tmass = 0;
  double ave_height;
  double height = kinetics->getMaxHeight();
  if (comm->me == 0) printf("max height = %e \n", height);
  bm1_output();
  // get biomass concentration (mol/L)
  for (int i = 0; i < nlocal; i++) {
    tmass += atom->rmass[i];
  }

  MPI_Allreduce(&tmass,&global_tmass,1,MPI_DOUBLE,MPI_SUM,world);
  //ave_height = cheight->compute_scalar();

  if (bm1cflag == 1) {
    if (global_tmass > 2e-12) {
      if (comm->me == 0) printf("tmass = %e \n\n", global_tmass);

      kinetics->monod->external_gflag = 0;
      kinetics->diffusion = this->diffusion;
      // solve mass balance in bulk liquid
      //kinetics->diffusion->bulkflag = 1;

      int k = 0;
      while(k < 10000) {
        kinetics->integration();
        for (int i = 1; i <= nnus; i++) {
          if (strcmp(bio->nuName[i], "o2") != 0) {
            diffusion->compute_bulk(i);
          }
        }
        bm1_output();
        k++;
      }
    }
  }

  // case 3, average biofilm thickness: Lf = 20 μm
  if (bm1cflag == 3) {
    if (global_tmass > 7.9e-14) {
      if (comm->me == 0) printf("tmass = %e \n\n", global_tmass);

      kinetics->monod->external_gflag = 0;
      kinetics->diffusion = this->diffusion;
      // solve mass balance in bulk liquid
      kinetics->diffusion->bulkflag = 1;

      int k = 0;
      while(k < 1000) {
        kinetics->integration();
        for (int i = 1; i <= nnus; i++) {
          if (strcmp(bio->nuName[i], "o2") != 0) {
            diffusion->compute_bulk(i);
          }
        }
        bm1_output();
        k++;
      }
    } else {
      kinetics->diffusion = NULL;
    }
  }
}

void FixVerify::bm1_output() {
  for (int nu = 1; nu <= nnus; nu++) {
    if (strcmp(bio->nuName[nu], "sub") == 0) {
      if (comm->me == 0) printf("S-sub-bulk = %e\n", kinetics->diffusion->nuBS[nu]);
      double s = get_ave_s_sub_base();
      if (comm->me == 0) printf("S-sub-base = %e\n", s);
    }

    if (strcmp(bio->nuName[nu], "o2") == 0) {
      if (comm->me == 0) printf("S-o2-bulk = %e\n", kinetics->diffusion->nuBS[nu]);
      double s = get_ave_s_o2_base();
      if (comm->me == 0) printf("S-o2-base = %e\n\n", s);
    }
  }
}

double FixVerify::get_ave_s_sub_base() {
  double ave_sub_s = 0;
  double global_ave_sub_s = 0;

  int nX = kinetics->subn[0] + 2;
  int nY = kinetics->subn[1] + 2;

  for (int grid = 0; grid < kinetics->diffusion->nXYZ; grid++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuName[nu], "sub") == 0) {
        int up = grid + nX * nY;

        if (kinetics->diffusion->xGrid[grid][2] < kinetics->zlo && !kinetics->diffusion->ghost[up]) {
          ave_sub_s += kinetics->diffusion->nuGrid[nu][grid];
        }
      }
    }
  }

  MPI_Allreduce(&ave_sub_s,&global_ave_sub_s,1,MPI_DOUBLE,MPI_SUM,world);

  return global_ave_sub_s/(kinetics->nx * kinetics->ny);
}

double FixVerify::get_ave_s_o2_base() {
  double ave_o2_s = 0;
  double global_ave_o2_s = 0;

  int nX = kinetics->subn[0] + 2;
  int nY = kinetics->subn[1] + 2;

  for (int grid = 0; grid < kinetics->diffusion->nXYZ; grid++) {
    for (int nu = 1; nu <= nnus; nu++) {
      if (strcmp(bio->nuName[nu], "o2") == 0) {
        int up = grid + nX * nY;

        if (kinetics->diffusion->xGrid[grid][2] < kinetics->zlo && !kinetics->diffusion->ghost[up]) {
          ave_o2_s += kinetics->diffusion->nuGrid[nu][grid];
        }
      }
    }
  }

  MPI_Allreduce(&ave_o2_s,&global_ave_o2_s,1,MPI_DOUBLE,MPI_SUM,world);

  return global_ave_o2_s/(kinetics->nx * kinetics->ny);
}

/* ---------------------------------------------------------------------- */

int FixVerify::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"demflag") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal fix_modify command");
    if (strcmp(arg[1],"no") == 0) demflag = 0;
    else if (strcmp(arg[1],"yes") == 0) demflag = 1;
    else error->all(FLERR,"Illegal demflag parameter (yes or no)");
    return 2;
  }
  return 0;
}

void FixVerify::remove_atom(double height) {
  AtomVecBio *avec = (AtomVecBio *) atom->style_match("bio");

  int i = 0;
  while (i < nlocal) {
    if (atom->x[i][2]+atom->radius[i] > height) {
      avec->copy(nlocal-1,i,1);
      nlocal--;
    } else {
      i++;
    }
  }
}

int FixVerify::reachHeight(double height) {
  const int nlocal = atom->nlocal;
  double * const * const x = atom->x;
  double * const r = atom->radius;
  int cout = 0;

  for (int i=0; i < nlocal; i++) {
    if((x[i][2] + r[i]) > height)
      cout ++;
  }

  int global_max;
  MPI_Allreduce(&cout, &global_max, 1, MPI_DOUBLE, MPI_MAX, world);

  return cout > 10 ? 1 : 0;
}
