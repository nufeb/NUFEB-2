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

#include "fix_bio_kinetics_energy.h"

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

FixKineticsEnergy::FixKineticsEnergy(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 4) error->all(FLERR,"Not enough arguments in fix kinetics/monod command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[3+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[3+i][2]);
  }

  kinetics = NULL;
  epsflag = 0;

}

/* ---------------------------------------------------------------------- */

FixKineticsEnergy::~FixKineticsEnergy()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixKineticsEnergy::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsEnergy::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
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
    } else if (strcmp(modify->fix[j]->style,"eps_extract") == 0) {
      epsflag = 1;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  EPSdens = input->variable->compute_equal(ivar[0]);

  bio = kinetics->bio;

  if (bio->nnus == 0)
    error->all(FLERR,"fix_kinetics/monod requires Nutrients input");
  else if (bio->catCoeff == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Catabolism Coeffs input");
  else if (bio->anabCoeff == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Anabolism Coeffs input");
  else if (bio->maintain == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Maintenance input");
  else if (bio->decay == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Decay input");
  else if (bio->ks == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Ks input");
  else if (bio->decayCoeff == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Decay Coeffs input");
  else if (bio->q == NULL)
    error->all(FLERR,"fix_kinetics/monod requires Consumption Rate input");

  nnus = bio->nnus;
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

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
  metabolism and atom update
------------------------------------------------------------------------- */
void FixKineticsEnergy::growth(double dt)
{
  mask = atom->mask;
  nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;
  type = atom->type;
  ntypes = atom->ntypes;

  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  maintain = bio->maintain;
  decay = bio->decay;

  radius = atom->radius;
  rmass = atom->rmass;
  outerMass = avec->outerMass;
  outerRadius = avec->outerRadius;

  nuS = kinetics->nuS;
  nuR = kinetics->nuR;
  DGRCat = kinetics->DRGCat;
  nuConv = kinetics->nuConv;
  gYield = kinetics->gYield;

  memory->create(gMonod,atom->ntypes+1,kinetics->bgrids,"kinetics/monod:gMonod");

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    #pragma omp parallel for
    for (int j = 0; j < kinetics->bgrids; j++) {
      gMonod[i][j] = -1;
      //minCatMonod[i][j] = -1;
    }
  }

  #pragma omp parallel for
  for (int i = 0; i < nall; i++) {
    //get new growth rate based on new nutrients
    double mass = biomass(i) * dt;
    //update bacteria mass, radius etc
    bio_update(mass, i);
  }

  memory->destroy(gMonod);
}

double FixKineticsEnergy::biomass(int i) {

  int t = type[i];
  int pos = position(i);

  if (bio->q[t] == 0 || DGRCat[t][pos] == 0) return 0;

  double qMet;       // specific substrate uptake rate for growth metabolism
  double qCat;       // specific substrate uptake rate for catabolism

  double bacMaint;    // specific substrate consumption required for maintenance
  double mu;          // specific biomass growth
  double mass;
  double invYield;

  double m = gMonod[t][pos];
  if (m < 0) {
    gMonod[t][pos] = grid_monod(pos, t, 1);
    qMet = bio->q[t] * gMonod[t][pos];
  } else {
    qMet = bio->q[t] * m;
  }

  qCat = qMet;
  bacMaint = maintain[t] / -DGRCat[t][pos];

  if (gYield[t][pos] != 0) invYield = 1/gYield[t][pos];
  else invYield = 0;

  for (int nu = 1; nu <= nnus; nu++) {
    if ((bio->nuType[nu] == 0 && bio->diffCoeff[nu] != 0)) {
      if (!nuConv[nu]) {
        double consume;

        if (1.2 * bacMaint < qCat) {
          double metCoeff = catCoeff[t][nu] * invYield + anabCoeff[t][nu];
          mu = gYield[t][pos] * (qMet - bacMaint);
          mass = mu * rmass[i];
          consume = mass  * metCoeff;
        } else if (qCat <= 1.2 * bacMaint && bacMaint <= qCat) {
          consume = catCoeff[t][nu] * gYield[t][pos] * bacMaint * rmass[i];
        } else {
          double f;
          if (bacMaint == 0) f = 0;
          else f = (bacMaint - qCat) / bacMaint;

          mass = -decay[t] * f * rmass[i];
          consume = -mass * bio->decayCoeff[t][nu] + catCoeff[t][nu] * gYield[t][pos] * qCat * rmass[i];
        }
        // convert biomass unit from kg to mol
        consume = consume * 1000 / 24.6;
        //calculate liquid consumption, mol/m3
        consume = consume / vol;

        nuR[nu][pos] += consume;
      }
    }
  }

  if (1.2 * bacMaint < qCat) {
    mu = gYield[t][pos] * (qMet - bacMaint);
    mass = mu * rmass[i];
  } else if (qCat <= 1.2 * bacMaint && bacMaint <= qCat) {
    mass = 0;
  } else {
    double f;
    if (bacMaint == 0) f = 0;
    else f = (bacMaint - qCat) / bacMaint;

    mass = -decay[t] * f * rmass[i];
  }

  return mass;
}

/* ----------------------------------------------------------------------
  get monod term w.r.t all nutrients
------------------------------------------------------------------------- */

double FixKineticsEnergy::grid_monod(int pos, int type, int ind)
{
  double monod = 1;

  //printf ("invYield = %e \n", invYield );
  for (int i = 1; i <= nnus; i++ ) {
    double s = nuS[i][pos];
    double ks = bio->ks[type][i];

    if (ks != 0) {
      if (s <= 0) return 0;
      monod *= s/(ks + s);
    }
  }

  return monod;
}

void FixKineticsEnergy::bio_update(double biomass, int i)
{
  double density;
  const double threeQuartersPI = (3.0/(4.0*MY_PI));
  const double fourThirdsPI = 4.0*MY_PI/3.0;
  const double third = 1.0/3.0;

  density = rmass[i] / (fourThirdsPI * radius[i] * radius[i] * radius[i]);
  rmass[i] = rmass[i] + biomass;

  //update mass and radius
  if (mask[i] == avec->maskHET) {
    //update HET radius
    radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
    //update mass and radius for EPS shell if PES production is on
    if (epsflag == 1) {
      outerMass[i] = fourThirdsPI * (outerRadius[i] * outerRadius[i] * outerRadius[i] - radius[i] * radius[i] * radius[i]) * EPSdens
          + biomass;
      outerRadius[i] = pow(threeQuartersPI * (rmass[i] / density + outerMass[i] / EPSdens), third);
    }
  } else if (mask[i] != avec->maskEPS && mask[i] != avec->maskDEAD){
    radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
    outerMass[i] = rmass[i];
    outerRadius[i] = radius[i];
  }
}


int FixKineticsEnergy::position(int i) {

  // get index of grid containing i
  int xpos = (atom->x[i][0] - xlo) / stepx + 1;
  int ypos = (atom->x[i][1] - ylo) / stepy + 1;
  int zpos = (atom->x[i][2] - zlo) / stepz + 1;
  int pos = (xpos - 1) + (ypos - 1) * nx + (zpos - 1) * (nx * ny);

  if (pos >= kinetics->bgrids) {
    printf("Too big! pos=%d   size = %i\n", pos, kinetics->bgrids);
  }

  return pos;
}
