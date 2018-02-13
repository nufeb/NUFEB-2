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

#include "fix_bio_fluid.h"

#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "force.h"
#include "pointers.h"
#include "update.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixFluid::FixFluid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix nufebFoam command");

  var = new char*[4];
  ivar = new int[4];

  for (int i = 0; i < 4; i++) {
    int n = strlen(&arg[3+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3+i][2]);
  }
}

FixFluid::~FixFluid()
{
  int i;
  for (i = 0; i < 4; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;
}

void FixFluid::init()
{

  for (int n = 0; n < 4; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix nufebFoam does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix nufebFoam is invalid style");
  }

  dem_steps = input->variable->compute_equal(ivar[0]);
  bio_steps = input->variable->compute_equal(ivar[1]);
  bio_dt = input->variable->compute_equal(ivar[2]);
  nloops = input->variable->compute_equal(ivar[3]);

  if (dem_steps < 0)
    error->warning(FLERR, "Negative input value for DEM steps, steps will be calculated and adjusted in nufebFoam");
  if (bio_steps < 0)
    error->all(FLERR, "Bio steps can not be negative");
  if (bio_dt < 0)
    error->all(FLERR, "Bio timestep can not be negative");
  if (nloops < 0)
    error->all(FLERR, "nufebFoam loop can not be negative");

  demflag = 0;
  dem_dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixFluid::setmask()
{
  return 0;
}
