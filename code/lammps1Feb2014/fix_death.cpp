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
  if (narg != 9) error->all(FLERR,"Illegal fix death command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix death command");

  var = new char*[4];
  ivar = new int[4];

  int i;
  for (i = 0; i < 4; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  seed = atoi(arg[8]);

  if (seed <= 0) error->all(FLERR,"Illegal fix death command: seed should be greater than 0");

  // Random number generator, same for all procs
  random = new RanPark(lmp,seed);
}

/* ---------------------------------------------------------------------- */

FixDeath::~FixDeath()
{
  delete random;
  int i;
  for (i = 0; i < 4; i++) {
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
  for (i = 0; i < 4; i++) {
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

  double R6 = input->variable->compute_equal(ivar[0]);
  double R7 = input->variable->compute_equal(ivar[1]);
  double R8 = input->variable->compute_equal(ivar[2]);
  double factor = input->variable->compute_equal(ivar[3]);

  double virtualMass = 0.0;
  double averageMass = 1e-13;

  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      double gHET = 0;
      double gAOB = 0;
      double gNOB = 0;
      if (mask[i] == 1) {
        gHET = 1;
      }
      if (mask[i] == 2) {
        gAOB = 1;
      }
      if (mask[i] == 3) {
        gNOB = 1;
      }
      virtualMass += ((gHET * R6) + (gAOB * R7) + (gNOB * R8)) * rmass[i];
  

    }

  }

  int kill = (random->uniform() * (nall - 1));

  int killed = 0;

  double VM = virtualMass/averageMass;

  fprintf(stdout, "Virtual Mass ratio: %f\n", VM);

  while (virtualMass > factor * averageMass) {
  	if (mask[kill] == groupbit) {
  		double gHET = 0;
        double gAOB = 0;
        double gNOB = 0;
        if (mask[kill] == 1) {
          gHET = 1;
        }
        if (mask[kill] == 2) {
          gAOB = 1;
        }
        if (mask[kill] == 3) {
          gNOB = 1;
        }
  		mask[kill] = 5;
  		virtualMass -= ((gHET * R6) + (gAOB * R7) + (gNOB * R8)) * rmass[kill];
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

  fprintf(stdout, "Killed: %i\n", killed);


  modify->addstep_compute(update->ntimestep + nevery);
}
