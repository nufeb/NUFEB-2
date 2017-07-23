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

#include "compute_bio_height.h"

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebHeight::ComputeNufebHeight(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal compute average height command");

  scalar_flag = 1;
  extscalar = 0;

  nx = atoi(arg[3]);
  ny = atoi(arg[4]);

  if (nx <= 0 || ny <= 0) error->all(FLERR,"Illegal compute nx or ny value");

  nxy = nx * ny;
}

/* ---------------------------------------------------------------------- */

ComputeNufebHeight::~ComputeNufebHeight()
{
   delete [] maxh;
}

void ComputeNufebHeight::init()
{
  maxh = new double[nxy]();

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
}

/* ---------------------------------------------------------------------- */
double ComputeNufebHeight::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **x = atom->x;

  scalar = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int pos = position(i);
      double z = x[i][2] + atom->radius[i];
      if (z > maxh[pos]) maxh[pos] = z;
    }
  }

  for (int i = 0; i < nxy; i++) {
      scalar += maxh[i] * stepx * stepy;
  }

  scalar = scalar/(xhi*yhi);
  return scalar;
}

/* ---------------------------------------------------------------------- */

int ComputeNufebHeight::position(int i) {

  // get index of grid containing i
  int xpos = (atom->x[i][0] - xlo) / stepx + 1;
  int ypos = (atom->x[i][1] - ylo) / stepy + 1;

  int pos = (xpos - 1) + (ypos - 1) * ny;

  if (pos >= nxy) {
     printf("Too big! pos=%d   size = %i\n", pos, nxy);
  }

  return pos;
}
