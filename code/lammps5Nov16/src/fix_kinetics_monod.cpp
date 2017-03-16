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

  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics/monod command");

  diffevery = force->inumeric(FLERR,arg[3]);
  if (diffevery < 0) error->all(FLERR,"Illegal fix kinetics/monod command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  kinetics = NULL;
  diffusion = NULL;

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
    delete [] minCatMonod[i];
  }

  delete [] gMonod;
  delete [] minCatMonod;
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
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  if (diffevery > 0 && diffusion == NULL)
    lmp->error->all(FLERR,"The fix diffusion command is required");

  EPSdens = input->variable->compute_equal(ivar[0]);

  bio = kinetics->bio;

  if (bio->nnus == 0)
    error->all(FLERR,"fix_kinetics requires # of Nutrients input");
  else if (bio->catCoeff == NULL)
    error->all(FLERR,"fix_kinetics requires Catabolism Coeffs input");
  else if (bio->anabCoeff == NULL)
    error->all(FLERR,"fix_kinetics requires Anabolism Coeffs input");
  else if (bio->nuChr == NULL)
    error->all(FLERR,"fix_kinetics requires Nutrient Charge input");
  else if (bio->maintain == NULL)
    error->all(FLERR,"fix_kinetics requires Maintenance input");
  else if (bio->decay == NULL)
    error->all(FLERR,"fix_kinetics requires Decay input");
  else if (bio->ks == NULL)
    error->all(FLERR,"fix_kinetics requires Ks input");
  else if (bio->decayCoeff == NULL)
    error->all(FLERR,"fix_kinetics requires Decay Coeffs input");

  ntypes = atom->ntypes;
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;
  gYield = kinetics->gYield;

  ngrids = nx * ny * nz;
  nnus = bio->nnus;
  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  maintain = bio->maintain;
  decay = bio->decay;

  gMonod = new double*[ntypes+1];
  minCatMonod = new double*[ntypes+1];

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    gMonod[i] = new double[ngrids];
    minCatMonod[i] = new double[ngrids];
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
  nuR = kinetics->nuR;
  DGRCat = kinetics->DRGCat;

  //initialization
  for (int i = 0; i <= ntypes; i++) {
    for (int j = 0; j < ngrids; j++) {
      gMonod[i][j] = -1;
      minCatMonod[i][j] = -1;
    }
  }

  if(diffevery != 0 && !(update->ntimestep % diffevery)) {
    //get nutrient consumption
    for (int i = 0; i < nall; i++) {
      //printf("mass = %e \n", rmass[i]);
      if (mask[i] & groupbit) {
        growth(i);
      }
    }
    //printf("Average R = %e \n", ar/ngrids);
    //solve diffusion
    diffusion->diffusion(0);
    //printf("R = %e \n", nuR[1][0]);

    for (int i = 0; i < nall; i++) {
      int pos = position(i);
      //get new growth rate based on new nutrients
      double biomass = growth(i) * update->dt;
      //update bacteria mass, radius etc
      bio_update(biomass, i);
    }

    //reset consumption
    for (int i = 0; i <= nnus; i++) {
      for (int j = 0; j < ngrids; j++) {
        nuR[i][j] = 0;
      }
    }
  } else {
    for (int i = 0; i < nall; i++) {
      if (mask[i] & groupbit) {
        double biomass = growth(i) * update->dt;
       // agr = agr + growthRate;
        //update bacteria mass, radius etc
        bio_update(biomass, i);
      }
    }
  }
  //if(!(update->ntimestep % 1000)) printf("Average growth rate = %e \n", agr/nall);
 // printf ("new mass[1] = %e \n", rmass[1]);
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

double FixKineticsMonod::growth(int i) {

  //calculate growth rate using minimum monod
  int t = type[i];
  int pos = position(i);
  double muMet;
  double muCat;

  double qmet = 0.0;
  double bacMaint;
  double mu;
  double rg;

  double m = gMonod[t][pos];
  if (m < 0) {
    gMonod[t][pos] = grid_monod(pos, t, 1);
    qmet = avec->atom_mu[i] * gMonod[t][pos];
  } else {
    qmet = avec->atom_mu[i] * m;
  }

  muMet = qmet;
// printf ("qmet = %e \n",qmet * 3600);
//
//  m = minCatMonod[t][pos];
//  if (m < 0) {
//    minCatMonod[t][pos] = grid_monod(pos, t, 2);
//    qmet = avec->atom_mu[i] * minCatMonod[t][pos];
//  } else {
//    qmet = avec->atom_mu[i] * m;
//  }

  muCat = qmet;
//  printf("muMet = %e \n", muCat * 3600);

  bacMaint = maintain[t] / -DGRCat[t][pos];
  //printf ("bacMaint = %e \n", bacMaint * 3600);
  mu = gYield[t][pos] * (muMet - bacMaint);
 // printf ("mu = %e \n", mu);

  for (int j = 1; j <= nnus; j++) {
    if (strcmp(bio->nuName[j], "h") != 0 && strcmp(bio->nuName[j], "h2o") != 0) {
      double consume;
      if (1.2 * bacMaint < muCat) {
        double invYield;
        if (gYield[t][pos] != 0)
          invYield = 1/gYield[t][pos];
        else
          invYield = 0;

        double metCoeff = catCoeff[t][j] * invYield + anabCoeff[t][j];
        //printf ("metCoeff %i = %e \n",j, metCoeff);
        rg = mu * rmass[i];
        consume = rg * metCoeff;
      } else if (muCat <= 1.2 * bacMaint && bacMaint <= muCat) {
        rg = 0;
        consume = catCoeff[t][j] * gYield[t][pos] * bacMaint * rmass[i];
      } else {
        double f = (bacMaint - muCat) / bacMaint;
        rg = -decay[t] * f * rmass[i];
        consume = -(rg) * bio->decayCoeff[t][j] + catCoeff[t][j] * gYield[t][pos] * muCat * rmass[i];
      }
      //printf ("consume = %e \n", 3600 * consume/24.6e-3);
      //calculate liquid consumption, convert from m3 to L, g to mol
      double uptake = consume / (vol * 24.6);
      //printf("uptake[%i] = %e vol = %e \n" , j, uptake, vol);
      nuR[j][pos] += uptake;
    }
  }
  //printf ("rg = %e \n", rg*3600);
  return rg;
}
//
///* ----------------------------------------------------------------------
//  get minimum monod term w.r.t all nutrients
//------------------------------------------------------------------------- */
//
//double FixKineticsMonod::minimal_monod(int pos, int type, int ind)
//{
//  vector<double> mon;
//  int size = 0;
//  double min = 0;
//  double invYield = 1/gYield[type][pos];
//  //printf ("invYield = %e \n", invYield );
//  for (int i = 1; i <= nnus; i++ ) {
//    if (strcmp(bio->nuName[i], "h") != 0 && strcmp(bio->nuName[i], "h2o") != 0) {
//      double coeff;
//      if (ind == 1) coeff = catCoeff[type][i] * invYield + anabCoeff[type][i];
//      else coeff = catCoeff[type][i];
//      //printf ("coeff = %e,ind = %i \n", coeff, ind);
//      if (coeff < 0) {
//        double v = nuS[i][pos]/(bio->ks[type][i] + nuS[i][pos]);
//        if (v != 0) {
//          mon.push_back(v);
//          size++;
//        }
//      }
//     }
//   }
//
//  if (mon.size() > 0)
//    min = *min_element(mon.begin(), mon.end());
//
//  //printf ("min = %e \n",nuS[1][pos]);
//  return min;
//}

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
  //printf ("growthRate = %e \n", growthRate);
  //printf ("old mass[i] = %e \n", rmass[i]);

  density = rmass[i] / (fourThirdsPI * radius[i] * radius[i] * radius[i]);
  rmass[i] = rmass[i] + biomass;
 // cout<<"before mass " <<growthRate << endl;
  //printf ("new mass[i] = %e \n", rmass[i]);
  //update mass and radius
  if (mask[i] == avec->maskHET) {
      outerMass[i] = fourThirdsPI * (outerRadius[i] * outerRadius[i] * outerRadius[i] - radius[i] * radius[i] * radius[i]) * EPSdens
      + biomass * nevery;

      outerRadius[i] = pow(threeQuartersPI * (rmass[i] / density + outerMass[i] / EPSdens), third);
      radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
  } else if (mask[i] != avec->maskEPS){
      radius[i] = pow(threeQuartersPI * (rmass[i] / density), third);
      outerMass[i] = 0.0;
      outerRadius[i] = radius[i];
  }
}
