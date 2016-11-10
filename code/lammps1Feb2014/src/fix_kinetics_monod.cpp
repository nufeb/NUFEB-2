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

#include <fix_kinetics_monod.h>
#include <fix_kinetics.h>
#include "atom_vec_bio.h"
#include <bio.h>
#include "math.h"
#include "string.h"
#include "stdlib.h"
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
#include "comm.h"
#include "domain.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsMonod::FixKineticsMonod(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 6) error->all(FLERR,"Not enough arguments in fix kinetics/monod command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/monod command");

  var = new char*[2];
  ivar = new int[2];

  for (int i = 0; i < 2; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
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

  for (int n = 0; n < 2; n++) {
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
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  bio = kinetics->bio;
  ntypes = atom->ntypes;

  rg = input->variable->compute_equal(ivar[0]);
  gvol = input->variable->compute_equal(ivar[1]);

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  ngrids = nx * ny * nz;
  nnus = bio->nnus;
  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  yield = bio->yield;

  metCoeff = kinetics->metCoeff;

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
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int* type = atom->type;

  nuS = kinetics->nuS;
  nuR = kinetics->nuR;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outerMass = avec->outerMass;
  double *outerRadius = avec->outerRadius;

  double** minMonod = new double*[ntypes+1];

  double density;
  const double threeQuartersPI = (3.0/(4.0*MY_PI));
  const double fourThirdsPI = 4.0*MY_PI/3.0;
  const double third = 1.0/3.0;

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    minMonod[i] = new double[ngrids];
    for (int j = 0; j < ngrids; j++) {
      minMonod[i][j] = -1;
    }
  }

  for (int i = 0; i <= nnus; i++) {
    for (int j = 0; j < ngrids; j++) {
      nuR[i][j] = 0;
    }
  }

  //Metabolism
  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      double growthRate = 0;
      double growthBac = 0;
      // get index of grid containing i
      int xpos = (atom->x[i][0] - xlo) / stepx + 1;
      int ypos = (atom->x[i][1] - ylo) / stepy + 1;
      int zpos = (atom->x[i][2] - zlo) / stepz + 1;
      int pos = (xpos - 1) + (ypos - 1) * ny + (zpos - 1) * (nx * ny);

      if (pos >= ngrids) {
         printf("Too big! pos=%d   size = %i\n", pos, ngrids);
      }

      //calculate growth rate using minimum monod
      int t = type[i];
      double monod = minMonod[t][pos];
      if (monod < 0) {
        minMonod[t][pos] = minimal_monod(pos, t);
        growthRate = avec->atom_growth[i] * minMonod[t][pos];
      } else {
        growthRate = avec->atom_growth[i] * monod;
      }
      //calculate amount of biomass formed
      growthBac = growthRate * atom->rmass[i];
     // cout << growthRate << endl;
      for (int i = 1; i <= nnus; i++) {
        double consume = metCoeff[t][i] * growthBac;

        if(bio->nuType[i] == 0) {
          //calculate liquid concentrations
          double sLiq = consume/vol*1000;
          //5.0000e-12
          nuR[i][pos] += sLiq;
        }
//          else if (atom->nuType[i] == 1) {
//          // calculate gas partial pressures
//          double pGas = consume * rg * temp / gvol;
//          nuR[i][pos] += pGas;
//        }
      }

      //update mass and radius
      if (mask[i] & groupbit) {
        double value = growthRate * update->dt;

        rmass[i] = rmass[i] * (1 + (value * nevery));
        density = rmass[i] / (fourThirdsPI * radius[i]*radius[i]*radius[i]);
        radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
      }
    }
  }

//  for (int i = 1; i <= ntypes; i++) {
//    cout << atom->typeName[i] << endl;
//      for (int j = 0; j < 5; j++) {
//          cout << atom->typeG[i][j] << " ";
//      }
//    cout << endl;
//  }

  for (int i = 0; i <= ntypes; i++) {
    delete [] minMonod[i];
  }
  delete [] minMonod;
}

/* ----------------------------------------------------------------------
  get minimum monod term w.r.t all nutrients
------------------------------------------------------------------------- */

double FixKineticsMonod::minimal_monod(int pos, int type)
{
  vector<double> mon;
  int size = 0;

  for (int i = 1; i <= nnus; i++ ) {
    if (matConsume[type][i] != 0 && bio->nuType[i] == 0) {
      double v = nuS[i][pos]/(bio->ks[type] + nuS[i][pos]);
      mon.push_back(v);
      size++;
    }
  }
  double min = *min_element(mon.begin(), mon.end());
  return min;
}
