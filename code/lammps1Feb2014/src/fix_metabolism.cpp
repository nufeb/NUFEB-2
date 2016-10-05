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

#include <fix_metabolism.h>
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
#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>
#include <iomanip>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;
using namespace Eigen;


/* ---------------------------------------------------------------------- */

FixMetabolism::FixMetabolism(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if (narg != 4) error->all(FLERR,"Not enough arguments in fix metabolism command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix metabolism command");

}

/* ---------------------------------------------------------------------- */

FixMetabolism::~FixMetabolism()
{
  memory->destroy(anabCoeff);
  memory->destroy(catCoeff);
}

/* ---------------------------------------------------------------------- */

int FixMetabolism::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetabolism::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  nnus = atom->nNutrients;
  ntypes = atom->ntypes;
  catCoeff = atom->catCoeff;
  anabCoeff = atom->anabCoeff;
  yield = atom->yield;
  vecConc = atom->vecConc;
  vecR = atom->vecR;

  metCoeff_calculus();

  nx = domain->nx;
  ny = domain->ny;
  nz = domain->nz;

  ngrids = nx * ny * nz;

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
//
//  //initialise cells
//  int grid = 0;
//  for (double i = xlo - (stepx/2); i < xhi + stepx; i += stepx) {
//    for (double j = ylo - (stepy/2); j < yhi + stepy; j += stepy) {
//      for (double k = zlo - (stepy/2); k < zhi + stepy; k += stepy) {
//        cellVol[cell] = xstep * ystep * zstep;
//      }
//    }
//  }
  //printf("stepX = %e \n", stepX);
}

/* ----------------------------------------------------------------------
  create metabolic coefficient for all microbial species
------------------------------------------------------------------------- */

void FixMetabolism::metCoeff_calculus()
{
  metCoeff = memory->create(metCoeff,ntypes+1,nnus+1,"atom:metCoeff");
  matCons = memory->create(matCons,ntypes+1,nnus+1,"atom:matCons");

  for (int i = 1; i < ntypes+1; i++) {
    for (int j = 1; j < nnus+1; j++) {
      if (strcmp(atom->nuName[j], "h") == 0 || strcmp(atom->nuName[j], "h2o") == 0) {
        metCoeff[i][j] = 0;
        matCons[i][j] = 0;
        continue;
      }
      double pThY = 1/yield[i];
      double coeff = catCoeff[i][j] * pThY + anabCoeff[i][j];
      metCoeff[i][j] = coeff;
      if (coeff < 0) {
        matCons[i][j] = 1;
      } else {
        matCons[i][j] = 0;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMetabolism::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  metabolism();
}

void FixMetabolism::metabolism()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int* type = atom->type;

  double** biomass = new double[ntypes+1];
  double** minMonod = new double*[ntypes+1];

  for (int i = 1; i <= ntypes; i++) {
    minMonod[i] = new double[ngrids];
    biomass[i] = new double[ngrids]();
    for (int j = 0; j < ngrids; j++) {
      minMonod[i][j] = -1;
    }
  }

  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      // get indices of grid containing atom
      int xpos = (atom->x[i][0] - xlo) / stepx + 1;
      int ypos = (atom->x[i][1] - ylo) / stepy + 1;
      int zpos = (atom->x[i][2] - zlo) / stepz + 1;
      int pos = (xpos - 1) + (ypos - 1) * ny + (zpos - 1) * (nx * ny);

      if (pos >= ngrids) {
         printf("Too big! pos=%d   size = %i\n", pos, ngrids);
      }

      //calculate growth rate
      int type = type[i];
      double monod = minMonod[type][pos];
      if (monod < 0) {
        monod[pos] = minimal_monod(pos, type);
      }
      double growth = atom->growth[type] * minMonod;

      //calculate amount of biomass formed
      biomass[type][pos] = biomass[type][pos] + growth*atom->mass[i];
    }
  }

  for (int i = 1; i <= nnus; i++) {
    //create consumption vector
    VectorXd r (ngrids);
    r.setZero();
    for (int j = 0; j < ngrids; j++) {
      //total consumption of all types per grid
      double sNu;
      for (int k = 1; k <= ntypes; k++) {
        //amount of nutrient consumed per type per grid
        double cons = metCoeff[k][i] * biomass[i][j];
        //calculate liquid concentrations
        if(atom->nuType == 0) {
          double sLiq = cons/vol;
          sNu += sLiq;
          // calculate gas partial pressures
        } else if (atom->nuType == 1) {
          double sGas = cons*uniGas*temp/vol;
          sNu += sGas;
        }
      }
      r(j) = sNu;
    }
    vecR[i] = r;
  }

  for (int i = 1; i <= ntypes; i++) {
    delete biomass[i];
    delete minMonod[i];
  }

  delete [] biomass;
  delete [] minMonod;
}

double FixMetabolism::minimal_monod(int pos, int type)
{
  double* mon = new double[];
  int size = 0;

  for (int i = 1; i <= nnus; i++ ) {
    if (matCons[type] != 0) {
      mon = new double [size+1];
      mon[size] = (vecConc[i][pos]) / (atom->ks[type] + vecConc[i][pos]);
      size++;
    }
  }
  double min = *min_element(mon, mon+size);
  delete [] mon;
  return min;
}
