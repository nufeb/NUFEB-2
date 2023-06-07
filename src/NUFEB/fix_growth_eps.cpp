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

#include "fix_growth_eps.h"

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

FixGrowthEPS::FixGrowthEPS(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/growth/eps command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/eps requires grid_style nufeb/chemostat");

  isub = -1;
  decay = 0.0;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/eps command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthEPS::update_cells()
{
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      reac[isub][i] += decay * dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthEPS::update_atoms()
{
  for (int i = 0; i < grid->ncells; i++) {
    grid->growth[igroup][i][0] = -decay;
  }

  update_atoms_coccus();
}
