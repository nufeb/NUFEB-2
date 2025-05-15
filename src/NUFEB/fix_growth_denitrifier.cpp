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

#include "fix_growth_denitrifier.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGrowthDenit::FixGrowthDenit(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 36)
    error->all(FLERR, "Illegal fix nufeb/growth/denit command. Expected at least 36 parameters");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/denit requires grid_style nufeb/chemostat");

  iss = -1;
  io2 = -1;
  ino3 = -1;
  ino2 = -1;
  ino = -1;
  in2o= -1;
  
  k_s2 = 0.0;
  k_s3 = 0.0;
  k_s4 = 0.0;
  k_s5 = 0.0;
  
  k_oh2 = 0.0;
  k_oh3 = 0.0;
  k_oh4 = 0.0;
  k_oh5 = 0.0;

  k_no3 = 0.0;
  k_no2 = 0.0;
  k_n2o = 0.0;
  k_no = 0.0;

  k_13no = 0.0;
  k_14no = 0.0;
  k_15no = 0.0;

  eta_g2 = 0.0;
  eta_g3 = 0.0;
  eta_g4 = 0.0;
  eta_g5 = 0.0;

  growth = 0.0;
  yield = 1.0;
  decay = 0.0;

  iss = grid->find(arg[3]);
  if (iss < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named S");
  k_s2 = utils::numeric(FLERR,arg[4],true,lmp);
  k_s3 = utils::numeric(FLERR,arg[5],true,lmp);
  k_s4 = utils::numeric(FLERR,arg[6],true,lmp);
  k_s5 = utils::numeric(FLERR,arg[7],true,lmp);

  io2 = grid->find(arg[8]);
  if (io2 < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named O2");
  k_oh2 = utils::numeric(FLERR,arg[9],true,lmp);
  k_oh3 = utils::numeric(FLERR,arg[10],true,lmp);
  k_oh4 = utils::numeric(FLERR,arg[11],true,lmp);
  k_oh5 = utils::numeric(FLERR,arg[12],true,lmp);

  ino3 = grid->find(arg[13]);
  if (ino3 < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named NO3");
  k_no3 = utils::numeric(FLERR,arg[14],true,lmp);
  
  ino2 = grid->find(arg[15]);
  if (ino2 < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named NO2");
  k_no2 = utils::numeric(FLERR,arg[16],true,lmp);
  
  in2o = grid->find(arg[17]);
  if (in2o < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named N2O");
  k_n2o = utils::numeric(FLERR,arg[18],true,lmp);

  ino = grid->find(arg[19]);
  if (ino < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named NO");
  k_no = utils::numeric(FLERR,arg[19],true,lmp);
  k_13no = utils::numeric(FLERR,arg[19],true,lmp);
  k_14no = utils::numeric(FLERR,arg[20],true,lmp);
  k_15no = utils::numeric(FLERR,arg[21],true,lmp);


  int iarg = 22;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g2") == 0) {
       eta_g2 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g3") == 0) {
       eta_g3 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g4") == 0) {
       eta_g4 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g5") == 0) {
       eta_g5 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/denit command. Did not recognize argument name. Expected either growth, yield, maintain, decay, or eta_g2, through eta_g5");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthDenit::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      //double tmp1 = growth * conc[ino2][i] / (no2_affinity + conc[ino2][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
      //double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

      //reac[ino2][i] -= 1 / yield * tmp1 * dens[igroup][i];
      //reac[io2][i] -= (1.15 - yield) / yield * tmp1 * dens[igroup][i] + tmp2 * dens[igroup][i];
      //reac[ino3][i] += 1 / yield * tmp1 * dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthDenit::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    //double tmp1 = growth * conc[ino2][i] / (no2_affinity + conc[ino2][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
    //double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

    //grid->growth[igroup][i][0] = tmp1 - tmp2 - decay;
  }

  update_atoms_coccus();
}
