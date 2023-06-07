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
#include "fix_density.h"
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "domain.h"
#include "group.h"
#include "update.h"
#include "atom_masks.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDensity::FixDensity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nufeb/density") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nufeb/density command");

  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

int FixDensity::setmask()
{
  int mask = 0;
  mask |= POST_PHYSICS_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixDensity::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "nevery") == 0) {
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery <= 0) error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixDensity::post_physics_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
}

/* ---------------------------------------------------------------------- */

void FixDensity::compute()
{
  for (int igroup = 0; igroup < group->ngroup; igroup++) {
    for (int i = 0; i < grid->ncells; i++)
      grid->dens[igroup][i] = 0.0;
  }

  double vol = grid->cell_size * grid->cell_size * grid->cell_size;
  // including ghost atoms because there can be atoms that moved inside the
  //   sub-domain and were not yet exchanged
  // forward communication garantees that we have the latest ghost positions
  //   which were updated during initial integrate
  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    if (atom->x[i][0] >= domain->sublo[0] && atom->x[i][0] < domain->subhi[0] &&
        atom->x[i][1] >= domain->sublo[1] && atom->x[i][1] < domain->subhi[1] &&
        atom->x[i][2] >= domain->sublo[2] && atom->x[i][2] < domain->subhi[2]) {
      int cell = grid->cell(atom->x[i]);

      double d = atom->biomass[i] * atom->rmass[i] / vol;
      grid->dens[0][cell] += d;
      for (int igroup = 0; igroup < group->ngroup; igroup++) {
        if (atom->mask[i] & group->bitmask[igroup])
          grid->dens[igroup][cell] += d;
      }
    }
  }
}
