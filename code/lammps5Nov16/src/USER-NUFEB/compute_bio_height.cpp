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

#include <stdlib.h> // atoi

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"
#include "neighbor.h"
#include "pair.h"
#include "force.h"

#include "compute_bio_height.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebHeight::ComputeNufebHeight(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg) {
  if (narg != 5)
    error->all(FLERR, "Illegal compute average height command");

  scalar_flag = 1;
  extscalar = 0;

  nx = atoi(arg[3]);
  ny = atoi(arg[4]);

  if (nx <= 0 || ny <= 0)
    error->all(FLERR, "Illegal compute nx or ny value");

  nxy = nx * ny;
}

/* ---------------------------------------------------------------------- */

ComputeNufebHeight::~ComputeNufebHeight() {
  delete[] maxh;
}

void ComputeNufebHeight::init() {
  stepx = domain->prd[0] / nx;
  stepy = domain->prd[1] / ny;

  // neighbor->cutneighmax is not yet initialized
  double cutneighmax = force->pair->cutforce + neighbor->skin;
  if (stepx > cutneighmax || stepy > cutneighmax)
    error->all(FLERR, "Grid step size for compute ave_height must be smaller than master list distance cutoff\n");

  grid = Grid<double, 2>(Box<double, 2>(domain->boxlo, domain->boxhi), { nx, ny });
  subgrid = Subgrid<double, 2>(grid, Box<double, 2>(domain->sublo, domain->subhi));
  ReduceGrid<ComputeNufebHeight>::setup();

  nxy = subgrid.cell_count();
  maxh = new double[nxy]();
}

/* ---------------------------------------------------------------------- */

double ComputeNufebHeight::compute_scalar() {
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double **x = atom->x;

  std::fill(maxh, maxh + nxy, 0);

  ReduceGrid<ComputeNufebHeight>::exchange();

  scalar = 0;
  if (is_bottom_most()) {
    for (int i = 0; i < nlocal + nghost; i++) {
      if ((mask[i] & groupbit) && subgrid.is_inside( { x[i][0], x[i][1] }) && x[i][2] >= domain->sublo[2] && x[i][2] < domain->subhi[2]) {
        int cell = subgrid.get_index( { x[i][0], x[i][1] });
        double z = x[i][2] + atom->radius[i] - domain->boxlo[2];
        maxh[cell] = MAX(maxh[cell], z);
      }
    }
    for (int i = 0; i < nxy; i++) {
      scalar += maxh[i] * stepx * stepy;
    }
  }

  scalar = scalar / (domain->prd[0] * domain->prd[1]);

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

  return scalar;
}

/* ---------------------------------------------------------------------- */

bool ComputeNufebHeight::is_bottom_most() const {
  return domain->sublo[2] == domain->boxlo[2];
}
