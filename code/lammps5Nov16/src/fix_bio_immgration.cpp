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

#include "fix_bio_immgration.h"

#include <string.h>
#include "random_park.h"
#include <math.h>
#include <stdlib.h>
#include <random>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"
#include "domain.h"

#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace boost;
using namespace MathConst;
using namespace math;
using namespace std;

/* ---------------------------------------------------------------------- */

FixImmigration::FixImmigration(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix shear command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix divide command: calling steps should be positive integer");

  seed = atoi(arg[5]);
  if (seed <= 0) error->all(FLERR,"Illegal fix divide command: seed should be greater than 0");

  // Random number generator, same for all procs
  random = new RanPark(lmp,seed);

  var = new char*[2];
  ivar = new int[2];

  for (int i = 0; i < 2; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }
}

FixImmigration::~FixImmigration()
{
  delete random;
  for (int i = 0; i < 2; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixImmigration::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixImmigration::init()
{
  for (int i = 0; i < 2; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix shear does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix shear is invalid style");
  }

  divMass = input->variable->compute_equal(ivar[0]);
  density = input->variable->compute_equal(ivar[1]);
}

/* ---------------------------------------------------------------------- */

void FixImmigration::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  //int r = round(random->uniform());
  //if (r)
  immgration();

}

/* ---------------------------------------------------------------------- */

void FixImmigration::immgration() {
  double **x = atom->x;
  double **rmass = atom->rmass;
  double *radius = atom->radius;

  // sample from metacommunity with uniform distribution
  mt19937 generator(seed);
  uniform_int_distribution<int> uniformIdis(1,atom->ntypes);

  int newTypes = uniformIdis(generator);

  //Randomise mass of immigrant (Gamma distribution)
  double newMass;
  double newRadius;
  double mMin = 0.75 * divMass;
  double mMax = 0.95 * divMass;
  double diff = mMax - mMin;

  double* xMass = new double[1000];

  for (int i = 0; i < 1000; i++) {
    xMass[i] = mMin + (diff / 1000) * (i+1);
  }
  double param[2];

  gamfit(xMass, 1000, param);
  gamma_distribution<double> gammadis(param[0],param[1]);

  newMass = gammadis(generator);
  newRadius = pow((3.0/(4.0*MY_PI)) * (newMass / density), 1.0/3.0);

  // specify random position for immigrant attachment
  double newX, newY, newZ;

  uniform_real_distribution<> uniformRdisX(xlo+newRadius*2, xhi);
  uniform_real_distribution<> uniformRdisY(ylo+newRadius*2, yhi);
  newX = uniformRdisX(generator);
  newY = uniformRdisY(generator);



  // Determine if position is overlapping existing cell


  delete []xMass;
}

/* ----------------------------------------------------------------------
  rewrite Matlab gamfit function in c++
------------------------------------------------------------------------- */

double* FixImmigration::gamfit(double *xx, int length, double *param) {
  double sum = 0;
  double sumlog = 0;
  double xbar;
  double s;
  double a, b;

  for (int i = 0; i < length; i++){
    sum += xx[i];
    sumlog += log(xx[i]);
  }

  xbar = sum / length;
  s = log(xbar)-sumlog / length;
  a = (3 - s + sqrt((s - 3)*(s - 3) + 24 * s))/(12 * s);

  for (int i = 1; i < 5; i++) {
    a = a - (log(a) - digamma(a) - s)/(1 / a - trigamma(a));
  }

  b = xbar / a;
  param[0] = a;
  param[1] = b;

  //printf("a = %e, b = %e \n", a, b);
  return param;
}


