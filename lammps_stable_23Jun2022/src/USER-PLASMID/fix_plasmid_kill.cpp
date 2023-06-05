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

#include "fix_plasmid_kill.h"

#include <string.h>
#include "atom.h"
#include "error.h"
#include "modify.h"
#include "group.h"
#include "update.h"
#include "atom_masks.h"
#include "fix_property_plasmid.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPlasmidKill::FixPlasmidKill(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal nufeb/plasmid/kill command");
  
  idead = group->find(arg[3]);
  if (idead < 0)
    error->all(FLERR, "Can't find group in fix nufeb/plasmid/kill");
  idead = 1 | group->bitmask[idead];

  fix_plasmid = nullptr;

  int ifix = modify->find_fix_by_style("^nufeb/property/plasmid");
  if (ifix < 0 ) error->all(FLERR,"Illegal nufeb/plasmid/kill command: "
      "requires fix nufeb/property/plasmid");
  fix_plasmid = (FixPropertyPlasmid *) modify->fix[ifix];
}

/* ---------------------------------------------------------------------- */

int FixPlasmidKill::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixPlasmidKill::modify_param(int narg, char **arg)
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

void FixPlasmidKill::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
}

/* ---------------------------------------------------------------------- */

void FixPlasmidKill::compute()
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
