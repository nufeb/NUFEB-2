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

#include <string.h>
#include "fix_death.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeath::FixDeath(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/death command");
  
  idead = group->find(arg[3]);
  if (idead < 0)
    error->all(FLERR, "Can't find group");
  diameter = force->numeric(FLERR, arg[4]);
}

/* ---------------------------------------------------------------------- */

int FixDeath::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeath::post_integrate()
{
  int *type = atom->type;
  int *mask = atom->mask;
  double *radius = atom->radius;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (radius[i] < 0.5 * diameter) {
        mask[i] = idead;
      }
    }
  }
}
