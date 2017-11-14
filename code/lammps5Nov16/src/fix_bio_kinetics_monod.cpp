/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include <math.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>

#include "atom.h"
#include "atom_vec_bio.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "math_const.h"
#include "memory.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "fix_bio_kinetics_monod.h"
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsMonod::FixKineticsMonod(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics/growth/monod command");

  var = new char*[2];
  ivar = new int[2];

  for (int i = 0; i < 2; i++) {
    int n = strlen(&arg[3+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[3+i][2]);
  }

  kinetics = NULL;
}

/* ---------------------------------------------------------------------- */

FixKineticsMonod::~FixKineticsMonod()
{
  int i;
  for (i = 0; i < 2; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixKineticsMonod::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;

  memory->destroy(species);
  memory->destroy(growrate);
}

/* ---------------------------------------------------------------------- */

void FixKineticsMonod::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 2; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics/monod does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics/monod is invalid style");
  }

  // register fix kinetics with this class
  kinetics = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  EPSdens = input->variable->compute_equal(ivar[0]);
  etaHET = input->variable->compute_equal(ivar[1]);

  bio = kinetics->bio;

  if (bio->nnus == 0)
    error->all(FLERR,"fix_kinetics/monod requires Nutrients input");
  else if (bio->maintain == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Maintenance input");
  else if (bio->decay == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Decay input");
  else if (bio->ks == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Ks input");
  else if (bio->yield == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Yield input");
  else if (bio->mu == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Growth Rate input");

  ntypes = atom->ntypes;
  nnus = bio->nnus;
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species =  memory->create(species,ntypes+1,"monod:species");
  growrate = memory->create(growrate,ntypes+1,2,kinetics->ngrids,"monod:growrate");

  //Get computational domain size
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

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  vol = stepx * stepy * stepz;

  init_param();
}

/* ----------------------------------------------------------------------
  initialize growth parameters
------------------------------------------------------------------------- */
void FixKineticsMonod::init_param()
{
  isub = io2 = inh4 = ino2 = ino3 = 0;
  ieps = idead = 0;

  // initialize nutrient
  for (int nu = 1; nu <= nnus; nu++) {
    if (strcmp(bio->nuName[nu], "sub") == 0) isub = nu;
    else if (strcmp(bio->nuName[nu], "o2") == 0) io2 = nu;
    else if (strcmp(bio->nuName[nu], "nh4") == 0) inh4 = nu;
    else if (strcmp(bio->nuName[nu], "no2") == 0) ino2 = nu;
    else if (strcmp(bio->nuName[nu], "no3") == 0) ino3 = nu;
    else error->all(FLERR,"unknow nutrient in fix_kinetics/kinetics/monod");

  }

  if (isub == 0) error->all(FLERR,"fix_kinetics/kinetics/monod requires nutrient substrate");
  if (io2 == 0) error->all(FLERR,"fix_kinetics/kinetics/monod requires nutrient o2");
  if (inh4 == 0) error->all(FLERR,"fix_kinetics/kinetics/monod requires nutrient nh4");
  if (ino2 == 0) error->all(FLERR,"fix_kinetics/kinetics/monod requires nutrient no2");
  if (ino3 == 0) error->all(FLERR,"fix_kinetics/kinetics/monod requires nutrient no3");

  // initialize type
  for (int i = 1; i <= atom->ntypes; i++) {
    if (strcmp(bio->typeName[i], "eps") == 0) {
      species[i] = 4;
      ieps = i;
    } else if (strcmp(bio->typeName[i], "dead") == 0) {
      species[i] = 5;
      idead = i;
    } else {
      // take first three char
      char *name = new char[4];
      strncpy(name, bio->typeName[i], 3);
      name[3] = 0;

      if (strcmp(name, "het") == 0) species[i] = 1;
      else if (strcmp(name, "aob") == 0) species[i] = 2;
      else if (strcmp(name, "nob") == 0) species[i] = 3;
      else error->all(FLERR,"unknow species in fix_kinetics/kinetics/monod");

      delete[] name;
    }
  }
  if (ieps == 0) (error->warning(FLERR,"EPS is not defined in fix_kinetics/kinetics/monod"));
}

/* ----------------------------------------------------------------------
  metabolism and atom update
------------------------------------------------------------------------- */
void FixKineticsMonod::growth(double dt)
{

  xtype = memory->create(xtype,ntypes+1, kinetics->bgrids,"monod:xtype");
  for (int i = 1; i <= ntypes; i++) {
    for (int grid = 0; grid < kinetics->bgrids; grid++){
      xtype[i][grid] = 0;
    }
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal;// + atom->nghost;
  int *type = atom->type;
  int ntypes = atom->ntypes;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outerMass = avec->outerMass;
  double *outerRadius = avec->outerRadius;

  double *mu = bio->mu;
  double *decay =  bio->decay;
  double *maintain = bio->maintain;
  double *yield = bio->yield;
  double **ks = bio->ks;

  double **nuS = kinetics->nuS;
  double **nuR = kinetics->nuR;

  bool *nuConv = kinetics->nuConv;
  double yieldEPS = 0;
  if (ieps != 0) yieldEPS = yield[ieps];


  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      int pos = kinetics->position(i);
      int t = type[i];
      double rmassCellVol = rmass[i] / vol;

      xtype[t][pos] += rmassCellVol;
    }
  }

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    for (int i = 1; i <= ntypes; i++) {
      int spec = species[i];

      // HET monod model
      if (spec == 1) {
        double R1 = mu[i] * (nuS[isub][grid] / (ks[i][isub] + nuS[isub][grid])) * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R4 = etaHET * mu[i] * (nuS[isub][grid] / (ks[i][isub] + nuS[isub][grid])) *
            (nuS[ino3][grid] / (ks[i][ino3] + nuS[ino3][grid])) * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R5 = etaHET * mu[i] * (nuS[isub][grid] / (ks[i][isub] + nuS[isub][grid])) *
            (nuS[ino2][grid] / (ks[i][ino2] + nuS[ino2][grid])) * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R6 = decay[i];

        double R10 = maintain[i] * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R13 = (1 / 2.86) * maintain[i] * etaHET * (nuS[ino3][grid] / (ks[i][ino3] + nuS[ino3][grid]))
            * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R14 = (1 / 1.71) * maintain[i] * etaHET * (nuS[ino2][grid] / (ks[i][ino2] + nuS[ino2][grid]))
            * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));

        if (!nuConv[isub]) nuR[isub][grid] +=  ((-1 / yield[i]) * ((R1 + R4 + R5) * xtype[i][grid]));
        //if (xtype[i][grid] != 0) printf("nuR = %e \n", xtype[i][grid]);
        if (!nuConv[io2]) nuR[io2][grid] +=  (-((1 - yield[i] - yieldEPS) / yield[i]) * R1 * xtype[i][grid]);
        if (!nuConv[ino2]) nuR[ino2][grid] +=  -(((1 - yield[i] - yieldEPS) / (1.17 * yield[i])) * R5 * xtype[i][grid]);
        if (!nuConv[ino3]) nuR[ino3][grid] +=  -(((1 - yield[i] - yieldEPS) / (2.86 * yield[i])) * R4 * xtype[i][grid]);
        if (!nuConv[io2]) nuR[io2][grid] += -(R10 * xtype[i][grid]);
        if (!nuConv[ino2]) nuR[ino2][grid] += -(R14 * xtype[i][grid]);
        if (!nuConv[ino3]) nuR[ino3][grid] += -(R13 * xtype[i][grid]);

        growrate[i][0][grid] = dt * (R1 + R4 + R5 - R6 - R10 - R13 - R14);
        growrate[i][1][grid] = dt * (yieldEPS / yield[i]) * (R1 + R4 + R5);
      } else if (spec == 2) {
        // AOB monod model
        double R2 = mu[i] * (nuS[inh4][grid] / (ks[i][inh4] + nuS[inh4][grid])) * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R7 = decay[i];
        double R11 = maintain[i] * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));

        if (!nuConv[io2]) nuR[io2][grid] +=  -(((3.42 - yield[i]) / yield[i]) * R2 * xtype[i][grid]);
        if (!nuConv[inh4]) nuR[inh4][grid] +=  -(1 / yield[i]) * R2 * xtype[i][grid];
        if (!nuConv[ino2]) nuR[ino2][grid] +=  (1 / yield[i]) * R2 * xtype[i][grid];
        if (!nuConv[io2]) nuR[io2][grid] += -(R11 * xtype[i][grid]);

        growrate[i][0][grid] = dt * (R2 - R7 - R11);
      } else if (spec == 3) {
        // NOB monod model
        double R3 = mu[i] * (nuS[ino2][grid] / (ks[i][ino2] + nuS[ino2][grid])) * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));
        double R8 = decay[i];
        double R12 = maintain[i] * (nuS[io2][grid] / (ks[i][io2] + nuS[io2][grid]));

        if (!nuConv[io2]) nuR[io2][grid] +=  - (((1.15 - yield[i]) / yield[i]) * R3 * xtype[i][grid]);
        if (!nuConv[ino2]) nuR[ino2][grid] +=  (1 / yield[i]) * R3 * xtype[i][grid];
        if (!nuConv[ino3]) nuR[ino3][grid] +=  -(1 / yield[i]) * R3 * xtype[i][grid];
        if (!nuConv[io2]) nuR[io2][grid] += -(R12 * xtype[i][grid]);

        growrate[i][0][grid] = dt * (R3 - R8 - R12);
      } else if (spec == 4) {
        // EPS monod model
        double R9 = decay[i];

        if (!nuConv[isub]) nuR[isub][grid] += R9 * xtype[i][grid];
        growrate[i][0][grid] = dt * -decay[i];
      } else if (spec == 5) {
        // DEAD monod model
        if (!nuConv[isub]) nuR[isub][grid] += (decay[i] * xtype[i][grid]);
        growrate[i][0][grid] = dt * -decay[i];
      }
    }
  }

  const double threeQuartersPI = (3.0 / (4.0 * MY_PI));
  const double fourThirdsPI = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      int pos = kinetics->position(i);

      double density = rmass[i] / (fourThirdsPI * radius[i] * radius[i] * radius[i]);
      rmass[i] = rmass[i] * (1 + growrate[t][0][pos]);

      if (species[t] == 1) {
        outerMass[i] = fourThirdsPI * (outerRadius[i] * outerRadius[i] * outerRadius[i]
                - radius[i] * radius[i] * radius[i]) * EPSdens + growrate[t][1][pos] * rmass[i];

        outerRadius[i] = pow(threeQuartersPI * (rmass[i] / density + outerMass[i] / EPSdens), third);
        radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
      } else {
        radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
        outerMass[i] = 0.0;
        outerRadius[i] = radius[i];
      }
    }
  }

  memory->destroy(xtype);
}
