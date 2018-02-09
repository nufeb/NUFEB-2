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

#include "error.h"
#include "force.h"
#include "pointers.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixFluid::FixFluid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 11) error->all(FLERR,"Illegal fix nufebFoam command");

  dem_steps = force->inumeric(FLERR,arg[4]);
  bio_steps = force->inumeric(FLERR,arg[6]);
  bio_dt = force->numeric(FLERR, arg[8]);
  nloops = force->inumeric(FLERR, arg[10]);
}

FixFluid::~FixFluid()
{
}

void FixFluid::init()
{
  demflag = 0;
  dem_dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixFluid::setmask()
{
  return 0;
}
