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

#include "fix_bio_swim.h"

#include <string.h>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixSwim::FixSwim(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix swim command");

  nevery = force->inumeric(FLERR,arg[3]);

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }
}

FixSwim::~FixSwim()
{
  for (int i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixSwim::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSwim::init()
{
  for (int i = 0; i < 1; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix shear does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix shear is invalid style");
  }

  rate = input->variable->compute_equal(ivar[0]);
}

/* ---------------------------------------------------------------------- */

void FixSwim::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  double **f = atom->f;
  double **x = atom->x;
  double **v = atom->v;
  double *radius = atom->radius;

  int c = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    int numneigh = 0;
    for(int j = 0; j < atom->nlocal; j++){
      if(i != j){
        double xd = atom->x[i][0] - atom->x[j][0];
        double yd = atom->x[i][1] - atom->x[j][1];
        double zd = atom->x[i][2] - atom->x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (atom->radius[i] + atom->radius[j] + 2e-7) * (atom->radius[i] + atom->radius[j] + 2e-7);

        if (rsq <= cut) {
          numneigh++;
        }
      }
    }
    if (numneigh == 0) {
      c++;
      f[i][2] += rate;
    }
  }

  //printf("time = %i , c = %i \n", update->ntimestep, c);
}
