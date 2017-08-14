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

#include "fix_bio_kinetics_thermo.h"

#include <math.h>
#include <string.h>
#include <cstdio>
#include <string>
#include <sstream>

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "memory.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "modify.h"
#include "pointers.h"
#include "variable.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsThermo::FixKineticsThermo(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Not enough arguments in fix kinetics/thermo command");

  fixY = 0;
  closeR = 0;

  if (strcmp(arg[3], "unfix") == 0) fixY = 1;
  else if (strcmp(arg[3], "fix") == 0) fixY = 0;
  else error->all(FLERR,"Illegal fix kinetics/thermo command: specify 'fix' or 'unfix' yield");

  if (strcmp(arg[4], "close") == 0) closeR = 1;
  else if (strcmp(arg[4], "open") == 0) closeR = 0;
  else error->all(FLERR,"Illegal fix kinetics/thermo command: specify 'open' or 'close' reactor");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[5+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[5+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsThermo::~FixKineticsThermo()
{
  memory->destroy(khV);
  memory->destroy(liq2Gas);
  memory->destroy(dG0);

  for (int i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::init()
{
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
    lmp->error->all(FLERR,"The fix kinetics command is required for kinetics/thermo styles");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics/thermo does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics/thermo is invalid style");
  }

  pressure = input->variable->compute_equal(ivar[0]);

  bio = kinetics->bio;

  if (bio->catCoeff == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires Catabolism Coeffs input");
  else if (bio->anabCoeff == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires Anabolism Coeffs input");
  else if (bio->nuGCoeff == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires Nutrient Energy input");
  else if (bio->typeGCoeff == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires Type Energy input");
  else if (fixY == 1 && bio->dissipation == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires Dissipation inputs for unfix yield");
  else if (closeR == 1 && bio->kLa == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires KLa input for closed reactor");
  else if (fixY == 1 && bio->eD == NULL)
    error->all(FLERR,"fix_kinetics/thermo requires eD input for unfix yield");

  temp = kinetics->temp;
  rth = kinetics->rth;

  nnus = bio->nnus;
  kLa = bio->kLa;
  nuGCoeff = bio->nuGCoeff;

  khV = memory->create(khV,nnus+1,"kinetics/thermo:khV");
  liq2Gas = memory->create(liq2Gas,nnus+1,"kinetics/thermo:liq2Gas");
  dG0 = memory->create(dG0,atom->ntypes+1,2,"kinetics/thermo:dG0");

  init_KhV();

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

  stepx = (xhi - xlo) / kinetics->nx;
  stepy = (yhi - ylo) / kinetics->ny;
  stepz = (zhi - zlo) / kinetics->nz;

  vol = stepx * stepy * stepz;
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_dG0()
{
  int *ngflag = bio->ngflag;
  int *tgflag = bio->tgflag;

  for (int i = 1; i <= atom->ntypes; i++) {
    dG0[i][0] = 0;
    dG0[i][1] = 0;

    for (int j = 1; j <= nnus; j++) {
      int flag = ngflag[j];
      dG0[i][0] += catCoeff[i][j] * nuGCoeff[j][flag];
      dG0[i][1] += anabCoeff[i][j] * nuGCoeff[j][flag];
    }

    int flag = tgflag[i];
    dG0[i][1] += typeGCoeff[i][flag];
//    printf("dG0[%i][0] = %e \n", i, dG0[i][0]);
//    printf("dG0[%i][1] = %e \n", i, dG0[i][1]);
  }
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_KhV()
{
  for (int i = 1; i < nnus + 1; i++) {
    khV[i] = 0;
    liq2Gas[i] = 0;
  }

  for (int i = 1; i < nnus + 1; i++) {
    char *nName;     // nutrient name
    nName = bio->nuName[i];

    if (bio->nuType[i] == 1) {
      char *lName = new char[strlen(nName)];     // corresponding liquid
      strncpy(lName, nName+1, strlen(nName));

      for (int j = 1; j < nnus + 1; j++) {
        char *nName2;     // nutrient name
        nName2 = bio->nuName[j];
        if (strcmp(lName, nName2) == 0) {
          liq2Gas[i] = j;
          if (nuGCoeff[j][0] > 10000) {
            khV[j] = exp((nuGCoeff[j][1] - nuGCoeff[i][1]) / (-rth * temp));
          } else {
            khV[j] = exp((nuGCoeff[j][0] - nuGCoeff[i][1]) / (-rth * temp));
          }
          break;
        }
      }
      delete [] lName;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixKineticsThermo::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
  thermodynamics
------------------------------------------------------------------------- */
void FixKineticsThermo::thermo()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *type = atom->type;
  ntypes = atom->ntypes;

  DRGCat = kinetics->DRGCat;
  DRGAn = kinetics->DRGAn;
  nuR = kinetics->nuR;
  nuS = kinetics->nuS;
  gYield = kinetics->gYield;
  activity = kinetics->activity;
  qGas = kinetics->qGas;

  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  typeGCoeff = bio->typeGCoeff;
  diss = bio->dissipation;

  dG0 = memory->grow(dG0,ntypes+1,2,"kinetics/thermo:dG0");

  init_dG0();

  double vRgT = kinetics->gVol * 1000 / (kinetics->gasTrans * kinetics->temp);
  double vg = vol * 1000;
  double rGas, rLiq;

  if (closeR == 0) {
    for (int j = 1; j <= nnus; j++) {
      if (bio->nuType[j] == 1) {
        kinetics->nuConv[j] = true;
      }
    }
  }

  for (int i = 0; i < kinetics->bgrids; i++) {
    // gas transfer
    if (closeR == 1){
      for (int j = 1; j <= nnus; j++) {
        nuR[j][i] = 0;
        qGas[j][i] = 0;
        if (bio->nuType[j] == 1) {
          //get corresponding liquid ID
          int liqID = liq2Gas[j];
          double gasT = 0;
          if (liqID != 0) {
            if (nuGCoeff[liqID][0] > 10000) {
              gasT = kLa[liqID] * (activity[liqID][1][i] / khV[liqID] - activity[j][0][i]);
            } else {
              gasT = kLa[liqID] * (activity[liqID][0][i] / khV[liqID] - activity[j][0][i]);
            }
            rGas = gasT;
            rLiq = -gasT * vRgT;
            // update nutrient consumption
            nuR[j][i] += rGas;
            nuR[liqID][i] += rLiq;
            //gas correction
            nuR[j][i] = nuR[j][i] * (nuS[j][i] * kinetics->gasTrans * temp / pressure + vg);
            //update gas concentration
            nuS[j][i] += nuR[j][i];
            qGas[j][i] = nuR[j][i] * kinetics->gasTrans * temp / pressure;
            if (nuS[j][i] < 0) nuS[j][i] = 0;
          }
        }
      }
    }

    for (int j = 1; j <= ntypes; j++) {
      double rthT = temp * rth;

      //Gibbs free energy of the reaction
      DRGCat[j][i] = dG0[j][0];  //catabolic energy values
      DRGAn[j][i] = dG0[j][1] + rthT;  //anabolic energy values

      for (int k = 1; k <= nnus; k++) {
        if (nuGCoeff[k][1] < 1e4) {
          double value = 0;
          int flag = bio->ngflag[k];
          if (activity[k][flag][i] == 0) value = 1e-20;
          else value = activity[k][flag][i];

          double dgr = rthT * log(value);
         // printf ("%e ",  value);
          DRGCat[j][i] += catCoeff[j][k] * dgr;
          DRGAn[j][i] += anabCoeff[j][k] * dgr;
         //if (k ==1 && j==7)printf("%e %s %e %s\n", DRGCat[j][i], bio->typeName[j], catCoeff[j][k], bio->nuName[k]);
        } else {
          error->all(FLERR,"nuGCoeff[1] is inf value");
        }
      }
     // printf("\n");
//      printf("DRGAn[%i][%i][0] = %e \n", i,j, DRGAn[j][i]);
//      printf("DRGCat[%i][%i][1] = %e \n", i,j, DRGCat[j][i]);

      //use catabolic and anabolic energy values to derive catabolic reaction equation
      if (fixY == 1){
        if (DRGCat[j][i] < 0) {
          gYield[j][i] = -(DRGAn[j][i] + diss[j]) / DRGCat[j][i] + bio->eD[j];
          if (gYield[j][i] != 0)
            gYield[j][i] = 1 / gYield[j][i];
        } else {
          gYield[j][i] = 0;
        }
      }
    }
  }
}
