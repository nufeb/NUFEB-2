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
#include "fix_death_diameter.h"
#include "atom.h"
#include "error.h"
#include "group.h"
#include "update.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeathDiameter::FixDeathDiameter(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/death/diameter command");
  
  idead = group->find(arg[3]);
  if (idead == -1)
    error->all(FLERR, "Can't find group in fix nufeb/death/diameter");
  idead = 1 | group->bitmask[idead];

  if (idead < 0)
    error->all(FLERR, "Can't find group");
  diameter = utils::numeric(FLERR,arg[4],true,lmp);

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "type") == 0) {
      tdead = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/death/diameter command");
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixDeathDiameter::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixDeathDiameter::modify_param(int narg, char **arg)
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

void FixDeathDiameter::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
}

/* ---------------------------------------------------------------------- */

void FixDeathDiameter::compute()
{
  int *mask = atom->mask;
  int *type = atom->type;
  double *radius = atom->radius;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (radius[i] < 0.5 * diameter) {
        mask[i] = idead;
        if (tdead > 0) type[i] = tdead;
      }
    }
  }
}
