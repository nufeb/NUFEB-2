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

#include "fix_property_cycletime.h"

#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyCycletime::FixPropertyCycletime(LAMMPS *lmp, int narg, char **arg) :
    FixProperty(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nufeb/property/cycletime command");

  create_attribute = 1;
  scalar_flag = 1;
  // use aprop if size_peratom_cols > 0
  size_peratom_cols = 2;
  grow_arrays(atom->nmax);

  seed = 0;
  max_time = 0;
  ran_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "seed") == 0) {
      seed = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      random = new RanPark(lmp, seed);
      ran_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "max_time") == 0 ) {
      max_time = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/property/cycletime command");
    }
  }

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (ran_flag > 0) {
        aprop[i][0] = random->uniform() * max_time;
        aprop[i][1] = aprop[i][0];
      } else {
        aprop[i][0] = 0.0;
        aprop[i][1] = 0.0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixPropertyCycletime::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ----------------------------------------------------------------------
   update cell age
------------------------------------------------------------------------- */

void FixPropertyCycletime::biology_nufeb()
{
  if (update->ntimestep % nevery) return;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      aprop[i][0] += update->dt;
      aprop[i][1] += update->dt;
    }
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyCycletime::set_arrays(int j)
{
  aprop[j][0] = 0.0;
}

/* ----------------------------------------------------------------------
   update array values of two daughter cells i, j
   called in fix_divide
------------------------------------------------------------------------- */
void FixPropertyCycletime::update_arrays(int i, int j)
{
  // when cell divides, record cell cycle time and reset to 0
  aprop[i][0] = 0.0;
  aprop[j][1] = aprop[i][1];
}

/* ----------------------------------------------------------------------
   compute average cell cycle time
------------------------------------------------------------------------- */
double FixPropertyCycletime::compute_scalar() {
  double result = 0.0;
  int n = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((atom->mask[i] & groupbit) && aprop[i][1] > 0.0) {
      result += aprop[i][1];
      n++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_INT, MPI_SUM, world);

  if (n > 0) result /= n;
  else result = 0;

  return result;
}


