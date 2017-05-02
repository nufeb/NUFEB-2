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

#include "fix_bio_kinetics_monod.h"

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
#include "fix_bio_diffusion.h"
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

FixKineticsMonod::FixKineticsMonod(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics/monod command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/monod command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  kinetics = NULL;
  diffusion = NULL;
  epsflag = 0;

}

/* ---------------------------------------------------------------------- */

FixKineticsMonod::~FixKineticsMonod()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  for (int i = 0; i < ntypes + 1; i++) {
    delete [] gMonod[i];
    //delete [] minCatMonod[i];
  }

  delete [] gMonod;
  //delete [] minCatMonod;
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

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics/monod does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics/monod is invalid style");
  }

  // register fix kinetics with this class
  kinetics = NULL;
  diffusion = NULL;
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style,"diffusion") == 0) {
      diffusion = static_cast<FixDiffusion *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style,"eps_extract") == 0) {
      epsflag = 1;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  if (nevery > 0 && diffusion == NULL)
    lmp->error->all(FLERR,"The fix diffusion command is required");

  EPSdens = input->variable->compute_equal(ivar[0]);

  bio = kinetics->bio;

  if (bio->nnus == 0)
    error->all(FLERR,"fix_kinetics/monod requires # of Nutrients input");
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

  ntypes = atom->ntypes;
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;
  gYield = kinetics->gYield;

  ngrids = nx * ny * nz;
  nnus = bio->nnus;
  maintain = bio->maintain;
  decay = bio->decay;

  gMonod = new double*[ntypes+1];
  //minCatMonod = new double*[ntypes+1];

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    gMonod[i] = new double[ngrids];
    //minCatMonod[i] = new double[ngrids];
  }

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
  DGRCat = kinetics->DRGCat;

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    for (int j = 0; j < ngrids; j++) {
      gMonod[i][j] = -1;
      //minCatMonod[i][j] = -1;
    }
  }

  if(nevery != 0 && !(update->ntimestep % nevery)) {
    //solve diffusion
    diffusion->diffusion();

    for (int i = 0; i < nall; i++) {
      int pos = position(i);
      //get new growth rate based on new nutrients
      double biomass = grow(i) * update->dt;
      //update bacteria mass, radius etc
      bio_update(biomass, i);
    }
  } else {
    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        double biomass = grow(i) * update->dt;
       // agr = agr + growthRate;
        //update bacteria mass, radius etc
        bio_update(biomass, i);
      }
    }
  }
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

double FixKineticsMonod::grow(int i) {

  int t = type[i];
  int pos = position(i);

  if (avec->atom_mu[i] == 0 || DGRCat[t][pos] == 0)
    return 0;

  double qMet;       // specific substrate uptake rate for growth metabolism
  double qCat;       // specific substrate uptake rate for catabolism

  double bacMaint;    // specific substrate consumption required for maintenance
  double mu;          // specific biomass growth
  double biomass;

  double m = gMonod[t][pos];
  if (m < 0) {
    gMonod[t][pos] = grid_monod(pos, t, 1);
    qMet = avec->atom_mu[i] * gMonod[t][pos];
  } else {
    qMet = avec->atom_mu[i] * m;
  }

  qCat = qMet;
  bacMaint = maintain[t] / -DGRCat[t][pos];

  for (int j = 1; j <= nnus; j++) {
    if (strcmp(bio->nuName[j], "h") != 0 && strcmp(bio->nuName[j], "h2o") != 0) {
      double consume;

      if (1.2 * bacMaint < qCat) {
        double invYield;
        if (gYield[t][pos] != 0)
          invYield = 1/gYield[t][pos];
        else
          invYield = 0;

        mu = gYield[t][pos] * (qMet - bacMaint);
        biomass = mu * rmass[i];
      } else if (qCat <= 1.2 * bacMaint && bacMaint <= qCat) {
        biomass = 0;
      } else {
        double f = (bacMaint - qCat) / bacMaint;
        biomass = -decay[t] * f * rmass[i];
      }
    }
  }
  //printf ("rg = %e \n", rg*3600);
  return biomass;
}

/* ----------------------------------------------------------------------
  get monod term w.r.t all nutrients
------------------------------------------------------------------------- */

double FixKineticsMonod::grid_monod(int pos, int type, int ind)
{
  double monod = 1;

  //printf ("invYield = %e \n", invYield );
  for (int i = 1; i <= nnus; i++ ) {
    double ks = bio->ks[type][i];
    double s = nuS[i][pos];

    if (ks != 0) {
      if (s == 1e-20) return 0;
      monod *= s/(ks + s);
    }
  }

  return monod;
}

void FixKineticsMonod::bio_update(double biomass, int i)
{
  double density;
  const double threeQuartersPI = (3.0/(4.0*MY_PI));
  const double fourThirdsPI = 4.0*MY_PI/3.0;
  const double third = 1.0/3.0;

  density = rmass[i] / (fourThirdsPI * radius[i] * radius[i] * radius[i]);
  rmass[i] = rmass[i] + biomass * nevery;

  //update mass and radius
  if (mask[i] == avec->maskHET) {
    //update HET radius
    radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
    //update mass and radius for EPS shell if PES production is on
    if (epsflag == 1) {
      outerMass[i] = fourThirdsPI * (outerRadius[i] * outerRadius[i] * outerRadius[i] - radius[i] * radius[i] * radius[i]) * EPSdens
      + biomass * nevery;
      outerRadius[i] = pow(threeQuartersPI * (rmass[i] / density + outerMass[i] / EPSdens), third);
    }
  } else if (mask[i] != avec->maskEPS && mask[i] != avec->maskDEAD){
    radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
    outerMass[i] = 0.0;
    outerRadius[i] = radius[i];
  }
}
