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

#include "fix_divide.h"

#include <string.h>
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixDivide::FixDivide(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  force_reneighbor = 1;
}


/* ---------------------------------------------------------------------- */

int FixDivide::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixDivide::modify_param(int narg, char **arg)
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

void FixDivide::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixDivide::post_neighbor()
{
  // reset reneighbour flag
  next_reneighbor = 0;
}
