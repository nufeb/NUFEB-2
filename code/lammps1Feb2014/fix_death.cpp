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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_death.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDeath::FixDeath(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix death command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix death command");

  var = new char*[2];
  ivar = new int[2];

  int i;
  for (i = 0; i < 2; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  seed = atoi(arg[6]);

  if (seed <= 0) error->all(FLERR,"Illegal fix death command: seed should be greater than 0");

  // Random number generator, same for all procs
  random = new RanPark(lmp,seed);
}

/* ---------------------------------------------------------------------- */

FixDeath::~FixDeath()
{
  delete random;
  int i;
  for (i = 0; i < 2; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixDeath::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDeath::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix death requires atom attribute diameter");

  int i;
  for (i = 0; i < 2; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix death does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix death is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void FixDeath::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  death();
}

/* ---------------------------------------------------------------------- */

void FixDeath::death()
{

  modify->clearstep_compute();

  double decay = input->variable->compute_equal(ivar[0]);
  double factor = input->variable->compute_equal(ivar[1]);

  double virtualMass = 0.0;
  double averageMass = 1e-16;

  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

  for (i = 0; i < nall; i++) {
  //  fprintf(stdout, "Id, Type, Mask: %i, %i, %i\n", i, type[i], mask[i]);
    if (mask[i] & groupbit) {
      virtualMass += (decay * rmass[i]);
    }

  }

  // fprintf(stdout, "Virtual Mass: %e\n", virtualMass);

  int kill = (random->uniform() * (nall - 1));

  int killed = 0;

  double VM = virtualMass/averageMass;

  //fprintf(stdout, "Virtual Mass ratio: %f\n", VM);

  while (virtualMass > (factor * averageMass)) {
    //fprintf(stdout, "Virtual Mass ratio: %f\n", VM);
  	if (mask[kill] & groupbit) {
      //fprintf(stdout, "Killed\n");
  		type[kill] = 5;
      mask[kill] = 33;
  		virtualMass -= (decay * rmass[kill]);
  		kill = (random->uniform() * (nall - 1));
  		killed ++;
  	}
  	else {
  		kill++;
  		if (kill == nall) {
  			kill = 0;
  		}
  	}
  }

  //fprintf(stdout, "Killed: %i\n", killed);


  modify->addstep_compute(update->ntimestep + nevery);
}
