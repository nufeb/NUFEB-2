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

#include "fix_mutate.h"

#include <string.h>
#include "atom.h"
#include "error.h"
#include "modify.h"
#include "group.h"
#include "update.h"
#include "atom_masks.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMutate::FixMutate(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg < 6)
    error->all(FLERR, "Illegal fix mutation/mutate command");

  imutant = group->find(arg[3]);
  if (imutant < 0)
    error->all(FLERR, "Can't find group in fix mutation/mutate");
  imutant = 1 | group->bitmask[imutant];

  prob = utils::numeric(FLERR,arg[4],true,lmp);
  if (prob > 1 || prob < 0)
    error->all(FLERR, "Illegal fix mutation/mutate command");

  seed = utils::inumeric(FLERR,arg[5],true,lmp);
  if (seed < 0)
    error->all(FLERR, "Illegal fix mutation/mutate command");

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
}

/* ---------------------------------------------------------------------- */

int FixMutate::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixMutate::modify_param(int narg, char **arg)
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

void FixMutate::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
}

/* ---------------------------------------------------------------------- */

void FixMutate::compute()
{
  int *mask = atom->mask;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (random->uniform() < prob) {
        mask[i] = imutant;
      }
    }
  }
}
