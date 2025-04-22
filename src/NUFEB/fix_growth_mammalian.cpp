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

#include <cstring>
#include "atom.h"
#include "error.h"

#include "fix_growth_mammalian.h"
#include "grid.h"
#include "group.h"
#include "modify.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGrowthMammalian::FixGrowthMammalian(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR, "Illegal fix nufeb/growth/mammalian command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/mammalian requires grid_style nufeb/chemostat");

  dynamic_group_allow = 1;

  isub = -1;
  ivirus = -1;

  sub_affinity = 0.0;

  growth = 0.0;
  yield = 1.0;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate name");
  sub_affinity = utils::numeric(FLERR,arg[4],true,lmp);
  if (sub_affinity <= 0)
    error->all(FLERR, "Substrate affinity must be greater than zero");

  ivirus = grid->find(arg[5]);
  if (ivirus < 0)
    error->all(FLERR, "Can't find virus name");

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/mammalian command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthMammalian::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      double tmp1 = growth * conc[isub][i] / (sub_affinity + conc[isub][i]);

      // nutrient utilization
      reac[isub][i] -= 1 / yield * tmp1 * dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthMammalian::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    double tmp1 = growth * conc[isub][i] / (sub_affinity + conc[isub][i]);

    grid->growth[igroup][i][0] = tmp1;
  }
}
