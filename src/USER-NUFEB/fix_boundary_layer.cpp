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
#include "update.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBoundaryLayer::FixBoundaryLayer(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix nufeb/boundary_layer command");

  xlo = xhi = ylo = yhi = zlo = zhi = 0;
  nlayers = 0;

  height = utils::numeric(FLERR,arg[3],true,lmp);
  if (height < 0) error->all(FLERR, "Illegal fix nufeb/boundary_layer command");

  nlayers = utils::inumeric(FLERR,arg[4],true,lmp);
  for (int i = 0; i < nlayers; i++) {
    if (strcmp(arg[5+i], "xlo") == 0)
      xlo = 1;
    else if (strcmp(arg[5+i], "xhi") == 0)
      xhi = 1;
    else if (strcmp(arg[5+i], "ylo") == 0)
      ylo = 1;
    else if (strcmp(arg[5+i], "yhi") == 0)
      yhi = 1;
    else if (strcmp(arg[5+i], "zlo") == 0)
      zlo = 1;
    else if (strcmp(arg[5+i], "zhi") == 0)
      zhi = 1;
    else
      error->all(FLERR, "Illegal fix nufeb/boundary_layer command");
  }
}

/* ---------------------------------------------------------------------- */

int FixBoundaryLayer::setmask()
{
  int mask = 0;
  mask |= REACTOR_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixBoundaryLayer::modify_param(int narg, char **arg)
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

void FixBoundaryLayer::reactor_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
}

/* ---------------------------------------------------------------------- */

void FixBoundaryLayer::compute()
{
  double maximum[3], minimum[3];

  compute_extremumx(maximum, minimum);

  for (int i = 0; i < 3; i++) {
    layerhi[i] = static_cast<int>((height + maximum[i]) / grid->cell_size) + 1;
    layerlo[i] = static_cast<int>((minimum[i] - height) / grid->cell_size);
    sublayerhi[i] = layerhi[i] - (grid->sublo[i] + 1);
    sublayerlo[i] = layerlo[i] - (grid->sublo[i] + 1);
  }
//  printf("box0 = %i box1=%i box2=%i \n",grid->subbox[0],grid->subbox[1],grid->subbox[2] );
//  printf("xlo = %i xhi = %i ylo = %i yhi = %i zlo = %i zhi = %i \n ", sublayerlo[0],sublayerhi[0],sublayerlo[1],sublayerhi[1],sublayerlo[2],sublayerhi[2]);
  // update bulk mask
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
        if (grid->mask[i] & X_NB_MASK)
          m |= X_NB_MASK;
        if (grid->mask[i] & X_PB_MASK)
          m |= X_PB_MASK;
        if (grid->mask[i] & Y_NB_MASK)
          m |= Y_NB_MASK;
        if (grid->mask[i] & Y_PB_MASK)
          m |= Y_PB_MASK;
        if (grid->mask[i] & Z_NB_MASK)
          m |= Z_NB_MASK;
        if (grid->mask[i] & Z_PB_MASK)
          m |= Z_PB_MASK;

        if ((xlo && x <= sublayerlo[0]) ||
            (xhi && x >= sublayerhi[0]) ||
            (ylo && y <= sublayerlo[1]) ||
            (yhi && y >= sublayerhi[1]) ||
            (zlo && z <= sublayerlo[2]) ||
            (zhi && z >= sublayerhi[2])) {
          m |= BLAYER_MASK;
        }

        if (!(m & GHOST_MASK) && !(m & BLAYER_MASK)) {
          m |= GRID_MASK;
        }

        mask[i] = m;
      }
    }
  }
}


/* ----------------------------------------------------------------------
 compute maximum and minimum atom position in the system w.r.t 6 surfaces
 ------------------------------------------------------------------------- */
void FixBoundaryLayer::compute_extremumx(double *maximumx, double *minimumx) {
  double maximumx_local[3], minimumx_local[3];
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;

  minimumx_local[0] = domain->prd[0];
  minimumx_local[1] = domain->prd[1];
  minimumx_local[2] = domain->prd[2];
  maximumx_local[0] = maximumx_local[1] = maximumx_local[2] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (xhi && x[i][0] > maximumx_local[0])
	    maximumx_local[0] = x[i][0];
      if (yhi && x[i][1] > maximumx_local[1])
	    maximumx_local[1] = x[i][1];
      if (zhi && x[i][2] > maximumx_local[2])
	    maximumx_local[2] = x[i][2];
      if (xlo && x[i][0] < minimumx_local[0])
	    minimumx_local[0] = x[i][0];
      if (ylo && x[i][1] < minimumx_local[1])
	    minimumx_local[1] = x[i][1];
      if (zlo && x[i][2] < minimumx_local[2])
    	minimumx_local[2] = x[i][2];
    }
  }

  MPI_Allreduce(&maximumx_local[0], &maximumx[0], 3, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&minimumx_local[0], &minimumx[0], 3, MPI_DOUBLE, MPI_MIN, world);
}

