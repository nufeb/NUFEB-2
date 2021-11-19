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
#include "fix_boundary_layer.h"
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBoundaryLayer::FixBoundaryLayer(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix nufeb/boundary_layer command");

  compute_flag = 1;

  boundary = nullptr;

  isub = grid->find(arg[3]);

  height = utils::numeric(FLERR,arg[4],true,lmp);
  if (height < 0) error->all(FLERR, "Illegal fix nufeb/boundary_layer command");
}

/* ---------------------------------------------------------------------- */

int FixBoundaryLayer::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixBoundaryLayer::modify_param(int narg, char **arg)
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
void FixBoundaryLayer::init()
{
  boundary = grid->boundary[isub];
}

/* ---------------------------------------------------------------------- */

void FixBoundaryLayer::post_integrate()
{
  if (compute_flag)
    compute();
}

/* ---------------------------------------------------------------------- */

void FixBoundaryLayer::compute()
{
  double maximum[3], minimum[3];

  compute_extremumx(maximum, minimum);

  for (int i = 0; i < 3; i++) {
    layerhi[i] = static_cast<int>((height + maximum[i]) / grid->cell_size) + 1;
    layerlo[i] = static_cast<int>((minimum[i] - height) / grid->cell_size) - 1;
    sublayerhi[i] = layerhi[i] - (grid->sublo[i] + 1);
    sublayerlo[i] = layerlo[i] - (grid->sublo[i] + 1);
  }

  // update mask
  int *mask = grid->mask;
  for (int z = 0; z < grid->subbox[2]; z++) {
    for (int y = 0; y < grid->subbox[1]; y++) {
      for (int x = 0; x < grid->subbox[0]; x++) {
	int i = x + y * grid->subbox[0] + z * grid->subbox[0] * grid->subbox[1];
	int m = 0;
	if (grid->mask[i] & GHOST_MASK)
	  m |= GHOST_MASK;
	if (grid->mask[i] & CORNER_MASK)
	  m |= CORNER_MASK;
	else {
	  if ((grid->sublo[0] < 0 && x == 0) ||
	      (boundary[0] == DIRICHLET && x <= sublayerlo[0]))
	    m |= X_NB_MASK;
	  if ((grid->subhi[0] > grid->box[0] && x == grid->subbox[0] - 1) ||
	      (boundary[1] == DIRICHLET && x >= sublayerhi[0]))
	    m |= X_PB_MASK;
	  if ((grid->sublo[1] < 0 && y == 0) ||
	      (boundary[2] == DIRICHLET && y <= sublayerlo[1]))
	    m |= Y_NB_MASK;
	  if ((grid->subhi[1] > grid->box[1] && y == grid->subbox[1] - 1) ||
	      (boundary[3] == DIRICHLET && y >= sublayerhi[1]))
	    m |= Y_PB_MASK;
	  if ((grid->sublo[2] < 0 && z == 0) ||
	      (boundary[4] == DIRICHLET && z <= sublayerlo[2]))
	    m |= Z_NB_MASK;
	  if ((grid->subhi[2] > grid->box[2] && z == grid->subbox[2] - 1) ||
	      (boundary[5] == DIRICHLET && z >= sublayerhi[2]))
	    m |= Z_PB_MASK;
	}
	mask[i] = m;
      }
    }
  }
}


/* ----------------------------------------------------------------------
 taking the plane as substratum, get the maximum and minimum height of atoms
 ------------------------------------------------------------------------- */
void FixBoundaryLayer::compute_extremumx(double *maximumx, double *minimumx) {
  double maximumx_local[3], minimumx_local[3];
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double *radius = atom->radius;

  minimumx_local[0] = domain->prd[0];
  minimumx_local[1] = domain->prd[1];
  minimumx_local[2] = domain->prd[2];
  maximumx_local[0] = maximumx_local[1] = maximumx_local[2] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (boundary[1] == DIRICHLET && x[i][0] > maximumx_local[0])
	maximumx_local[0] = x[i][0];
      if (boundary[3] == DIRICHLET && x[i][1] > maximumx_local[1])
	maximumx_local[1] = x[i][1];
      if (boundary[5] == DIRICHLET && x[i][2] > maximumx_local[2])
	maximumx_local[2] = x[i][2];
      if (boundary[0] == DIRICHLET && x[i][0] < minimumx_local[0])
	minimumx_local[0] = x[i][0];
      if (boundary[2] == DIRICHLET && x[i][1] < minimumx_local[1])
	minimumx_local[1] = x[i][1] - radius[i];
      if (boundary[4] == DIRICHLET && x[i][2] < minimumx_local[2])
	minimumx_local[2] = x[i][2];
    }
  }

  MPI_Allreduce(&maximumx_local[0], &maximumx[0], 3, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&minimumx_local[0], &minimumx[0], 3, MPI_DOUBLE, MPI_MIN, world);
}

