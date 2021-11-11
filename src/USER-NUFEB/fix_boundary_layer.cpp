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
#include <cmath>
#include "fix_boundary_layer.h"
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "domain.h"
#include "group.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBoundaryLayer::FixBoundaryLayer(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix nufeb/boundary_layer command");

  xyl_flag = xyh_flag = yzl_flag = yzh_flag = xzl_flag = xzh_flag = 0;

  compute_flag = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "xyl") == 0) {
      xyl_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "xyh") == 0) {
      xyh_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "xzl") == 0) {
      xzl_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "xzh") == 0) {
      xzh_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "yzl") == 0) {
      yzl_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg], "yzh") == 0) {
      yzh_flag = 1;
      iarg++;
    } else {
      error->all(FLERR, "Illegal fix nufeb/boundary_layer command");
    }
  }

  nlayers = iarg-3;
}

/* ---------------------------------------------------------------------- */

int FixBoundaryLayer::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixBoundaryLayer::modify_param(int narg, char **arg)
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
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixBoundaryLayer::post_integrate()
{
  if (compute_flag)
    compute();
}

/* ---------------------------------------------------------------------- */

void FixBoundaryLayer::compute()
{

}
