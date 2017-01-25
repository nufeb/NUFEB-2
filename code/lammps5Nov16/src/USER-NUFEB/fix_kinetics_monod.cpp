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

#include "fix_kinetics_monod.h"

#include <math.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>

#include "atom.h"
#include "atom_vec_bio.h"
#include "bio.h"
#include "domain.h"
#include "error.h"
#include "fix_kinetics.h"
#include "fix_diffusion.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "math_const.h"
#include "memory.h"
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

  if (narg != 8) error->all(FLERR,"Not enough arguments in fix kinetics/monod command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/monod command");
  diffevery = force->inumeric(FLERR,arg[4]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/monod command");

  var = new char*[3];
  ivar = new int[3];

  for (int i = 0; i < 3; i++) {
    int n = strlen(&arg[5+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[5+i][2]);
  }

  kinetics = NULL;
  diffusion = NULL;
}

/* ---------------------------------------------------------------------- */

FixKineticsMonod::~FixKineticsMonod()
{
  int i;
  for (i = 0; i < 3; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  for (int i = 0; i <= ntypes; i++) {
    delete [] minMonod[i];
  }
  delete [] minMonod;

  memory->destroy(matConsume);
}

/* ---------------------------------------------------------------------- */

int FixKineticsMonod::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsMonod::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 3; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics/monod does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics/monod is invalid style");
  }

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style,"diffusion") == 0) {
      diffusion = static_cast<FixDiffusion *>(lmp->modify->fix[j]);
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  if (diffevery > 0 && diffusion == NULL)
    lmp->error->all(FLERR,"The fix diffusion command is required");

  bio = kinetics->bio;
  ntypes = atom->ntypes;

  rg = input->variable->compute_equal(ivar[0]);
  gvol = input->variable->compute_equal(ivar[1]);
  EPSdens = input->variable->compute_equal(ivar[2]);

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  ngrids = nx * ny * nz;
  nnus = bio->nnus;
  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  yield = bio->yield;

  metCoeff = kinetics->metCoeff;
  minMonod = new double*[ntypes+1];

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    minMonod[i] = new double[ngrids];
  }

  create_metaMatrix();

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
}

/* ----------------------------------------------------------------------
  calculate metabolic coefficient for all microbial species
------------------------------------------------------------------------- */

void FixKineticsMonod::create_metaMatrix()
{
  matConsume = memory->create(matConsume,ntypes+1,nnus+1,"atom:matCons");

  for (int i = 1; i <= ntypes; i++) {
    for (int j = 1; j <= nnus; j++) {
     // cout << "type = " << atom->nuType[j] << endl;
      if (strcmp(bio->nuName[j], "h") == 0 || strcmp(bio->nuName[j], "h2o") == 0) {
        metCoeff[i][j] = 0;
        matConsume[i][j] = 0;
        continue;
      }

      double pThY = 1/yield[i];
      double coeff = catCoeff[i][j] * pThY + anabCoeff[i][j];

      metCoeff[i][j] = coeff;
      if (coeff < 0) {
        matConsume[i][j] = 1;
      } else {
        matConsume[i][j] = 0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixKineticsMonod::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  monod();
}

/* ----------------------------------------------------------------------
  metabolism and atom update
------------------------------------------------------------------------- */
void FixKineticsMonod::monod()
{
  mask = atom->mask;
  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;
  type = atom->type;

  radius = atom->radius;
  rmass = atom->rmass;
  outerMass = avec->outerMass;
  outerRadius = avec->outerRadius;

  nuS = kinetics->nuS;
  nuR = kinetics->nuR;
 // double agr = 0;
 // double ar = 0;

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    for (int j = 0; j < ngrids; j++) {
      minMonod[i][j] = -1;
    }
  }

  if(diffevery != 0 && !(update->ntimestep % diffevery)) {
    //initialize consumption
    for (int i = 0; i <= nnus; i++) {
      for (int j = 0; j < ngrids; j++) {
        nuR[i][j] = 0;
      }
    }

    //get nutrient consumption
    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        int t = type[i];
        int pos = position(i);
        double growthRate = growth_rate(i, pos);
        //update bacteria mass, radius etc
        for (int j = 1; j <= nnus; j++) {
          double consume = metCoeff[t][j] * growthRate * rmass[i];
          if(bio->nuType[j] == 0) {
            //calculate liquid concentrations, convert from m3 to L, kg to mol
            double uptake = consume / (vol * 1000 * 2.46e-2);
            nuR[j][pos] += uptake;
           // ar += uptake;
          }
        }
      }
    }
    //printf("Average R = %e \n", ar/ngrids);
    //solve diffusion
    diffusion->diffusion(0);

    for (int i = 0; i < nall; i++) {
      int pos = position(i);
      //get new growth rate based on new nutrients
      double growthRate = growth_rate(i, pos) * update->dt;
      //agr = agr + growthRate;
      //update bacteria mass, radius etc
      bio_update(growthRate, i);
    }
  } else {
    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        int pos = position(i);
        double growthRate = growth_rate(i, pos) * update->dt;
       // agr = agr + growthRate;
        //update bacteria mass, radius etc
        bio_update(growthRate, i);
      }
    }
  }
  //if(!(update->ntimestep % 1000)) printf("Average growth rate = %e \n", agr/nall);
}

int FixKineticsMonod::position(int i) {

  // get index of grid containing i
  int xpos = (atom->x[i][0] - xlo) / stepx + 1;
  int ypos = (atom->x[i][1] - ylo) / stepy + 1;
  int zpos = (atom->x[i][2] - zlo) / stepz + 1;
  int pos = (xpos - 1) + (ypos - 1) * ny + (zpos - 1) * (nx * ny);

  if (pos >= ngrids) {
     printf("Too big! pos=%d   size = %i\n", pos, ngrids);
  }

  return pos;
}

double FixKineticsMonod::growth_rate(int i, int pos) {

  //calculate growth rate using minimum monod
  int t = type[i];
  double growthRate = 0.0;

  double monod = minMonod[t][pos];
  if (monod < 0) {
    minMonod[t][pos] = minimal_monod(pos, t);
    growthRate = avec->atom_mu[i] * minMonod[t][pos];
  } else {
    growthRate = avec->atom_mu[i] * monod;
  }

  return growthRate;
}

/* ----------------------------------------------------------------------
  get minimum monod term w.r.t all nutrients
------------------------------------------------------------------------- */

double FixKineticsMonod::minimal_monod(int pos, int type)
{
  vector<double> mon;
  int size = 0;
  double min = 0;

  for (int i = 1; i <= nnus; i++ ) {
    if (matConsume[type][i] != 0 && bio->nuType[i] == 0) {
      double v = nuS[i][pos]/(bio->ks[type] + nuS[i][pos]);
      mon.push_back(v);
      size++;
    }
  }

  if (mon.size() > 0)
    min = *min_element(mon.begin(), mon.end());

  return min;
}

void FixKineticsMonod::bio_update(double growthRate, int i)
{
  double density;
  const double threeQuartersPI = (3.0/(4.0*MY_PI));
  const double fourThirdsPI = 4.0*MY_PI/3.0;
  const double third = 1.0/3.0;

  density = rmass[i] / (fourThirdsPI * radius[i] * radius[i] * radius[i]);
  rmass[i] = rmass[i] * (1 + (growthRate * nevery));
 // cout<<"before mass " <<growthRate << endl;

  //update mass and radius
  if (mask[i] == avec->maskHET) {
      outerMass[i] = fourThirdsPI * (outerRadius[i] * outerRadius[i] * outerRadius[i] - radius[i] * radius[i] * radius[i]) * EPSdens
      + growthRate * nevery * rmass[i];

      outerRadius[i] = pow(threeQuartersPI * (rmass[i] / density + outerMass[i] / EPSdens), third);
      radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
  }
  else if (mask[i] != avec->maskEPS){
      radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
      outerMass[i] = 0.0;
      outerRadius[i] = radius[i];
  }
}
