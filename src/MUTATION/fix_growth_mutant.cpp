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

#include "fix_growth_mutant.h"

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

FixGrowthMutant::FixGrowthMutant(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 6)
    error->all(FLERR, "Illegal fix mutation/growth/mutant command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix mutation/growth/mutant requires grid_style nufeb/chemostat");

  isub = -1;
  iinh = -1;

  sub_affinity = 0.0;

  growth = 0.0;
  ic50 = 1.0;
  nl = 1.0;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate name");
  sub_affinity = utils::numeric(FLERR,arg[4],true,lmp);

  iinh = grid->find(arg[5]);
  if (iinh < 0)
    error->all(FLERR, "Can't find inhibitor name");

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "gamma") == 0) {
      gamma = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ic50") == 0) {
      ic50 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "nl") == 0) {
      nl = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }  else {
      error->all(FLERR, "Illegal fix mutation/growth/mutant command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthMutant::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      double tmp1 = growth * (conc[isub][i] / (sub_affinity + conc[isub][i])) / (1 + pow((conc[iinh][i] / ic50), nl));

      reac[isub][i] += -gamma * tmp1 * dens[igroup][i];
      reac[iinh][i] += 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthMutant::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    double tmp1 = growth * (conc[isub][i] / (sub_affinity + conc[isub][i])) / (1 + pow ((conc[iinh][i] / ic50),nl));

    grid->growth[igroup][i][0] = tmp1;
  }

  update_atoms_coccus();
}
