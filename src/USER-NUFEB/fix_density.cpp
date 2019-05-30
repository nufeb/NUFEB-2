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

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDensity::FixDensity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"nufeb/density") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix nufeb/density command");

  dynamic_group_allow = 1;
  compute_flag = 1;
}

/* ---------------------------------------------------------------------- */

int FixDensity::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixDensity::modify_param(int narg, char **arg)
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

void FixDensity::post_integrate()
{
  if (!compute_flag)
    return;
  
  for (int igroup = 0; igroup < group->ngroup; igroup++) {
    for (int i = 0; i < grid->ncells; i++)
      grid->dens[igroup][i] = 0.0;
  }

  double lo[3];
  double hi[3];
  for (int i = 0; i < 3; i++) {
    lo[i] = (grid->sublo[i] + 1) * grid->cell_size + domain->boxlo[i];
    hi[i] = (grid->subhi[i] - 1) * grid->cell_size + domain->boxlo[i];
  }

  double vol = grid->cell_size * grid->cell_size * grid->cell_size;
  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    if (atom->x[i][0] >= lo[0] && atom->x[i][0] < hi[0] &&
	atom->x[i][1] >= lo[1] && atom->x[i][1] < hi[1] &&
	atom->x[i][2] >= lo[2] && atom->x[i][2] < hi[2]) {
      int cell = grid->cell(atom->x[i]);
      double d = atom->rmass[i] / vol;
      grid->dens[0][cell] += d;
      for (int igroup = 0; igroup < group->ngroup; igroup++)
	if (atom->mask[i] & group->bitmask[igroup])
	  grid->dens[igroup][cell] += d;
    }
  }
}
