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

#include <cstdio>
#include <cstring>
#include "fix_monod.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMonod::FixMonod(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  compute_flag = 1;
  reaction_flag = 1;
  growth_flag = 1;
  dt = 1.0;
}

/* ---------------------------------------------------------------------- */

int FixMonod::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "compute") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	compute_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	compute_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "reaction") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	reaction_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	reaction_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "growth") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) {
	growth_flag = 1;
      } else if (strcmp(arg[iarg+1], "no") == 0) {
	growth_flag = 0;
      } else {
	error->all(FLERR, "Illegal fix_modify command");
      }
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixMonod::init()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixMonod::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixMonod::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMonod::post_integrate()
{
  if (compute_flag)
    compute();
}
