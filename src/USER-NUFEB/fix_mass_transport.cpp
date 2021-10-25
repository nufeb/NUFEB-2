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

#include "fix_mass_transport.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "error.h"
#include "comm.h"

#include "grid.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMassTransport::FixMassTransport(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix nufeb/gas_liquid command");

  compute_flag = 1;

  isub = -1;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate name");
}

/* ---------------------------------------------------------------------- */

int FixMassTransport::modify_param(int narg, char **arg)
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

int FixMassTransport::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMassTransport::post_integrate()
{
  if (compute_flag)
    compute();
}

/* ---------------------------------------------------------------------- */

void FixMassTransport::compute()
{
  double **conc = grid->conc;
  double **reac = grid->reac;

  // compute average substrate consumption
  double ave_reac = 0.0;
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      ave_reac += grid->reac[isub][i];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &ave_reac, 1, MPI_DOUBLE, MPI_SUM, world);
  ave_reac /= (grid->box[0] * grid->box[1] * grid->box[2]);

  // update substrate concentration
  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      grid->conc[isub][i] -= ave_reac;
    }
  }
}


