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

#include <cstdio>
#include <cstring>
#include <cmath>

#include "fix_biomass_balance.h"
#include "error.h"
#include "grid.h"
#include "grid_masks.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBiomassBalance::FixBiomassBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!grid->reactor_flag)
    error->all(FLERR,"Fix nufeb/biomass_balanc requires nufeb/reactor grid style");
  if (narg < 10)
    error->all(FLERR,"Illegal fix nufeb/biomass_balance command");

  q = 0.0;
  rvol = 1.0;
  af = 0.0;
  nsubs = 0;
  xy = 1.0;

  isub = NULL;
  inlet = NULL;

  if (strcmp(arg[3], "q") != 0)
    error->all(FLERR,"Illegal fix nufeb/biomass_balance command: q");
  q = force->numeric(FLERR, arg[4]);

  if (strcmp(arg[5], "rvol") != 0)
    error->all(FLERR,"Illegal fix nufeb/biomass_balance command: rvol");
  rvol = force->numeric(FLERR, arg[6]);
  if (rvol <= 0)
    error->all(FLERR, "Bioreactor volume cannot be less or equal to 0");

  if (strcmp(arg[7], "af") != 0)
    error->all(FLERR,"Illegal fix nufeb/biomass_balance command: af");
  af = force->numeric(FLERR, arg[8]);
  if (af < 0)
    error->all(FLERR, "Biofilm surface area cannot be less than 0");

  if (strcmp(arg[9], "xy") != 0)
    error->all(FLERR,"Illegal fix nufeb/biomass_balance command: xy");
  xy = force->numeric(FLERR, arg[10]);
  if (xy <= 0)
    error->all(FLERR, "Domain surface area cannot be less or equal to 0");

  int iarg = 11;
  nsubs = narg - iarg;

  isub = memory->create(isub, nsubs, "nufeb/biomass_balance:isub");
  for (int i = 0; i < nsubs; i++) {
    while (iarg < narg) {
      int ind = grid->find(arg[iarg]);
      if (ind < 0)
	error->all(FLERR, "Can't find substrate name");
      isub[i] = ind;
      iarg++;
    }
  }

  compute_flag = 1;
}

/* ---------------------------------------------------------------------- */

FixBiomassBalance::~FixBiomassBalance()
{

  if (copymode) return;
  memory->destroy(isub);
  memory->destroy(inlet);
}

/* ---------------------------------------------------------------------- */

void FixBiomassBalance::init()
{
  inlet = memory->create(inlet, grid->nsubs, "nufeb/biomass_balance:inlet");

  for (int i = 0; i < grid->nsubs; i++) {
    inlet[i] = grid->bulk[i];
  }
}

/* ---------------------------------------------------------------------- */

int FixBiomassBalance::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixBiomassBalance::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compute") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	compute_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	compute_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixBiomassBalance::post_integrate()
{
  if (compute_flag)
    compute();
}

/* ---------------------------------------------------------------------- */

void FixBiomassBalance::compute()
{
  double **reac = grid->reac;
  double *bulk = grid->bulk;
  double tot_r, global_tot_r;
  int ind;
  double vol = grid->cell_size * grid->cell_size * grid->cell_size;

  for (int i = 0; i < nsubs; i++) {
    tot_r = 0;
    global_tot_r = 0;
    ind = isub[i];

    for (int j = 0; j < grid->ncells; j++) {
      if (!(grid->mask[j] & GHOST_MASK)) {
	tot_r += reac[ind][j];
      }
    }
    MPI_Allreduce(&tot_r, &global_tot_r, 1, MPI_DOUBLE, MPI_SUM, world);
    // solve for the biomass balance in bulk liquid
    bulk[ind] += ((q / rvol) * (inlet[ind] - bulk[ind]) + ((af * global_tot_r * vol) / (rvol * xy))) * update->dt;
  }
}
