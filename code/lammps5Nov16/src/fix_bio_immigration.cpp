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

#include "fix_bio_immigration.h"

#include <string.h>
#include "random_park.h"
#include <math.h>
#include <stdlib.h>
#include <random>
#include <stdio.h>

#include "atom_vec_bio.h"
#include "fix_bio_kinetics.h"
#include "bio.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"
#include "domain.h"
#include "modify.h"

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
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 8) error->all(FLERR,"Illegal fix immigration command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix immigration command: calling steps should be positive integer");
  if (strcmp(arg[6], "on") == 0) zflag = 1;
  else if (strcmp(arg[6], "off") == 0) zflag = 0;
  seed = atoi(arg[7]);
  if (seed <= 0) error->all(FLERR,"Illegal fix immigration command: seed should be greater than 0");

  // Random number generator, same for all procs
  random = new RanPark(lmp,seed);
  kinetics = NULL;

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

  bio = avec->bio;

  for (int i = 0; i < atom->nlocal; i++)   atom->omega[i][0] = atom->type[i];

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
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
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixImmigration::init()
{
  for (int i = 0; i < 2; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix immigration does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix immigration is invalid style");
  }

  // register fix kinetics with this class
  kinetics = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  divMass = input->variable->compute_equal(ivar[0]);
  density = input->variable->compute_equal(ivar[1]);
}

/* ---------------------------------------------------------------------- */

void FixImmigration::post_integrate()
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  int r = round(random->uniform());
  if (r) immgration();
}

/* ---------------------------------------------------------------------- */

void FixImmigration::immgration() {
  // sample from metacommunity with uniform distribution
  int sampleT = 0;
  mt19937 generator(update->ntimestep);

  if (atom->ntypes == 1 && avec->typeEPS == 1) return;
  else if (atom->ntypes == 1 && avec->typeDEAD == 1) return;
  else if (atom->ntypes == 2 && avec->typeDEAD != 0  && avec->typeEPS != 0) return;

  // get sample type id
  do {
  uniform_int_distribution<int> uniformIdis(1,atom->ntypes);
  sampleT = uniformIdis(generator);
  } while (sampleT == avec->typeDEAD || sampleT == avec->typeEPS);

  // create new type
  char *oldName = bio->typeName[sampleT];
  char *newName = create_type_name(oldName);
  bio->create_type(newName);
  // new type id
  int newT = atom->ntypes;

  // take type parameters from sample type
  if (bio->yield != NULL) bio->yield[newT] = bio->yield[sampleT];
  if (bio->maintain != NULL) bio->maintain[newT] = bio->maintain[sampleT];
  if (bio->decay != NULL) bio->decay[newT] = bio->decay[sampleT];
  if (bio->dissipation != NULL) bio->dissipation[newT] = bio->dissipation[sampleT];
  if (bio->anabCoeff != NULL) memcpy(bio->anabCoeff[newT], bio->anabCoeff[sampleT], (bio->nnus + 1) * sizeof(double));
  if (bio->catCoeff != NULL) memcpy(bio->catCoeff[newT], bio->catCoeff[sampleT], (bio->nnus + 1) * sizeof(double));
  if (bio->decayCoeff != NULL) memcpy(bio->decayCoeff[newT], bio->decayCoeff[sampleT], (bio->nnus + 1) * sizeof(double));
  if (bio->typeGCoeff != NULL) memcpy(bio->typeGCoeff[newT], bio->typeGCoeff[sampleT], 5 * sizeof(double));
  if (bio->tgflag != NULL) bio->tgflag[newT] = bio->tgflag[sampleT];
  if (bio->typeChr != NULL) bio->typeChr[newT] = bio->typeChr[sampleT];
  if (bio->eD != NULL) bio->eD[newT] = bio->eD[sampleT];
  kinetics->grow();
//  if (atom->ntypes > 5)
//  test(5, 2);
  // take mu from gamma distribution
  // parameter estimation
  double muMin, muMax;
  double newMu;
  double scf = 0.05;
  double muDiff;

  double *xMu = new double[1000];

  muMin = bio->mu[sampleT] - scf * bio->mu[sampleT];
  muMax = bio->mu[sampleT] + scf * bio->mu[sampleT];
  muDiff = muMax - muMin;

  for (int i = 0; i < 1000; i++) {
    xMu[i] = muMin + (muDiff / 1000) * (i+1);
  }

  double muParam[2];
  gamfit(xMu, 1000, muParam);
  gamma_distribution<double> gammaMu(muParam[0],muParam[1]);
  newMu = gammaMu(generator);
  if (bio->mu != NULL) bio->mu[newT] = newMu;

  // take ks from gamma distribution
  for (int nu = 1; nu < bio->nnus+1; nu++) {
    if (bio->ks[sampleT][nu] != 0) {
      double ksMin, ksMax;
      double newKs;
      double scf2 = 0.05;
      double ksDiff;
      double *xKs = new double[1000];

      ksMin = bio->ks[sampleT][nu] - scf2 * bio->ks[sampleT][nu];
      ksMax = bio->ks[sampleT][nu] + scf2 * bio->ks[sampleT][nu];
      ksDiff = ksMax - ksMin;

      for (int i = 0; i < 1000; i++) {
        xKs[i] = ksMin + (ksDiff / 1000) * (i+1);
      }

      double ksParam[2];
      gamfit(xKs, 1000, ksParam);
      gamma_distribution<double> gammaKs(ksParam[0],ksParam[1]);
      newKs = gammaKs(generator);

      if (bio->ks != NULL) {
        if (bio->ks > 0)  bio->ks[newT][nu] = newKs;
        else  bio->ks[newT][nu] = 0;
      }

      delete []xKs;

    } else bio->ks[newT][nu] = 0;
  }

  //Randomise mass of immigrant (Gamma distribution)
  double newMass;
  double newRadius;
  double mMin = 0.75 * divMass;
  double mMax = 0.95 * divMass;
  double mDiff = mMax - mMin;

  double *xMass = new double[1000];

  for (int i = 0; i < 1000; i++) {
    xMass[i] = mMin + (mDiff / 1000) * (i+1);
  }
  double mParam[2];

  gamfit(xMass, 1000, mParam);
  gamma_distribution<double> gammaMass(mParam[0],mParam[1]);

  newMass = gammaMass(generator);
  newRadius = pow((3.0/(4.0*MY_PI)) * (newMass / density), 1.0/3.0);

  // specify random position for immigrant attachment
  double newX, newY, newZ;
  double* coord = new double[3];

  uniform_real_distribution<> uniformRdisX(xlo+newRadius*2, xhi-newRadius*2);
  uniform_real_distribution<> uniformRdisY(ylo+newRadius*2, yhi-newRadius*2);
  newX = uniformRdisX(generator);
  newY = uniformRdisY(generator);
  if (zflag == 0)
    newZ = find_z(newX, newY, newRadius);
  else
    newZ = zhi - newRadius;

  coord[0] = newX;
  coord[1] = newY;
  coord[2] = newZ;

  // sample atom
  int i = 0;
  for (i = 0; i < atom->nlocal; i++) {
    if (atom->type[i] == sampleT) break;
  }

  find_maxid();
  atom->avec->create_atom(newT,coord);
  // new atom
  int n = atom->nlocal-1;

  atom->tag[n] = maxtag_all+1;
  atom->mask[n] = atom->mask[i];
  atom->image[n] = atom->image[i];

  atom->v[n][0] = atom->v[i][0];
  atom->v[n][1] = atom->v[i][1];
  atom->v[n][2] = -10;

  atom->f[n][0] = atom->f[i][0];
  atom->f[n][1] = atom->f[i][1];
  atom->f[n][2] = -10;

  atom->omega[n][0] = atom->omega[i][0];
  atom->omega[n][1] = atom->omega[i][1];
  atom->omega[n][2] = atom->omega[i][2];

  atom->torque[n][0] = atom->torque[i][0];
  atom->torque[n][1] = atom->torque[i][1];
  atom->torque[n][2] = atom->torque[i][2];

  atom->rmass[n] = newMass;
  avec->outerMass[n] = newMass;

  atom->radius[n] = newRadius;
  avec->outerRadius[n] = newRadius;

  avec->atom_mu[n] = newMu;

  atom->natoms++;

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;

//  for (int i = 1; i<atom->ntypes+1; i++) {
//    printf("type = %s   ", bio->typeName[i]);
//  }
//  printf("\n");
//
//  for (int i = 0; i < atom->nlocal; i++) {
//    printf("i=%i, x=%e, y=%e, z=%e, typename = %s, type = %i \n", i, atom->x[i][0], atom->x[i][1], atom->x[i][2], bio->typeName[atom->type[i]], atom->type[i]);
//  }

  delete []xMass;
  delete []xMu;
  delete []coord;
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

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

char* FixImmigration::create_type_name(char *oldName){
  //printf("oldName = %s \n", oldName);
  int size = strlen(oldName);
  int no = 0;
  // take last integer
  char *p = oldName;
  while (*p) { // While there are more characters to process...
      if (isdigit(*p)) no = strtol(p, &p, 10); // Read a number, ...
      else p++;
  }

  int digit;
  if (no == 0) digit = 0;
  else digit= floor (log10 (abs (no))) + 1;
  int strsize = size - digit;
  char *str = new char[strsize+1];
  strncpy(str, oldName, strsize);
  str[strsize] = '\0';
  // maximum name id
  int max = 0;

  for (int i = 1; i < atom->ntypes+1; i++) {
    if (strstr(bio->typeName[i], str)) {
      int noi = 0;
      char *p = bio->typeName[i];
      while (*p) { // While there are more characters to process...
          if (isdigit(*p)) noi = strtol(p, &p, 10); // Read a number, ...
          else p++;
      }

      if (noi > max) max = noi;
    }
  }
  max++;

  int newDigit = floor (log10 (abs (max))) + 1;
  char *newNo = new char[newDigit+1];
  sprintf(newNo, "%i", max);
  char *newName = new char[strsize + newDigit + 1];
  strcpy(newName, str);
  strcat(newName, newNo);
  //printf("newName2 = %s \n", newName);

  delete []str;
  delete []newNo;

  return newName;
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */
double FixImmigration::find_z(double x1, double y1, double r){
  double z = r;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    double r = atom->radius[i];
    double** x = atom->x;

    if ((x[i][0] + r < x1) && (x[i][0] - r > x1) && (x[i][1] + r < y1) && (x[i][1] - r > y1))
      if (x[i][2] + r > z)
        z = x[i][2] + r;
  }

  return z;
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixImmigration::find_maxid()
{
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
}


/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixImmigration::test(int newT, int sampleT)
{
  printf("newT = %i, sampleT = %i \n", newT, sampleT);
  printf("newT = %s, sampleT = %s \n", bio->typeName[newT], bio->typeName[sampleT]);
  printf("nYield = %e, sYield = %e \n", bio->yield[newT], bio->yield[sampleT]);
  printf("newCat \n");
  for (int i = 1; i<=bio->nnus; i++) {
    printf("%e  ", bio->catCoeff[newT][i]);
  }
  printf("\n");
  printf("sampleCat \n");
  for (int i = 1; i<=bio->nnus; i++) {
    printf("%e  ", bio->catCoeff[sampleT][i]);
  }
  printf("\n");
}


