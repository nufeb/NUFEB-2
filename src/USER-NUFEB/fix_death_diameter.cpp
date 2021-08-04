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
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeathDiameter::FixDeathDiameter(LAMMPS *lmp, int narg, char **arg) :
  FixDeath(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/death command");

  compute_flag = 1;
  
  idead = group->find(arg[3]);
  if (idead < 0)
    error->all(FLERR, "Can't find group");
  diameter = utils::numeric(FLERR,arg[4],true,lmp);
}

/* ---------------------------------------------------------------------- */

void FixDeathDiameter::compute()
{
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
