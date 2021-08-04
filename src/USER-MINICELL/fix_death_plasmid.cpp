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
#include "atom.h"
#include "error.h"
#include "modify.h"
#include "group.h"
#include "atom_masks.h"
#include "fix_death_plasmid.h"
#include "fix_property_plasmid.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDeathPlasmid::FixDeathPlasmid(LAMMPS *lmp, int narg, char **arg) :
  FixDeath(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal nufeb/death/plasmid command");

  compute_flag = 1;
  
  idead = group->find(arg[3]);
  if (idead < 0)
    error->all(FLERR, "Can't find group");

  fix_plasmid = nullptr;;

  int ifix = modify->find_fix_by_style("^nufeb/property/plasmid");
  if (ifix < 0 ) error->all(FLERR,"Illegal nufeb/death/plasmid command: "
      "requires fix nufeb/property/plasmid");
  fix_plasmid = (FixPropertyPlasmid *) modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */

void FixDeathPlasmid::compute()
{
  int *mask = atom->mask;
  double *radius = atom->radius;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (!(int)fix_plasmid->vprop[i]) {
        mask[i] = idead;
      }
    }
  }
}
