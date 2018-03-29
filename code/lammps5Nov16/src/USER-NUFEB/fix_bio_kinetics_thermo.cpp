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

FixKineticsThermo::FixKineticsThermo(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  if (narg < 3)
    error->all(FLERR, "Not enough arguments in fix kinetics/thermo command");

  // default values
  yflag = 0;
  rflag = 0;
  pressure = 1;

  int iarg = 3;
  while (iarg < narg){
    if (strcmp(arg[iarg],"pressure") == 0) {
      pressure = force->numeric(FLERR, arg[iarg + 1]);
      if (pressure < 0.0)
        error->all(FLERR, "Illegal fix kinetics/thermo command: pressure");
      iarg += 2;
    } else if (strcmp(arg[iarg],"yield") == 0) {
      if (strcmp(arg[iarg+1],"fix") == 0) yflag = 0;
      else if (strcmp(arg[iarg+1],"unfix") == 0) yflag = 1;
      else error->all(FLERR,"Illegal yield parameter:fix or unfix");
      iarg += 2;
    } else if (strcmp(arg[iarg],"reactor") == 0) {
      if (strcmp(arg[iarg+1],"close") == 0) rflag = 0;
      else if (strcmp(arg[iarg+1],"open") == 0) rflag = 1;
      else error->all(FLERR,"Illegal yield parameter:open or close");
    } else
      error->all(FLERR, "Illegal fix kinetics/thermo command");
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsThermo::~FixKineticsThermo() {
  memory->destroy(khv);
  memory->destroy(liqtogas);
  memory->destroy(dgzero);
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::init() {
  // register fix kinetics with this class
  kinetics = NULL;
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for kinetics/thermo styles");

  bio = kinetics->bio;

  if (bio->catCoeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Catabolism Coeffs input");
  else if (bio->anabCoeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Anabolism Coeffs input");
  else if (bio->nuGCoeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Nutrient Energy input");
  else if (bio->typeGCoeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Type Energy input");
  else if (yflag == 1 && bio->dissipation == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Dissipation inputs for unfix yield");
  else if (rflag == 1 && bio->kLa == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires KLa input for closed reactor");
  else if (yflag == 1 && bio->eD == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires eD input for unfix yield");

  int nnus = bio->nnus;

  khv = memory->create(khv, nnus + 1, "kinetics/thermo:khv");
  liqtogas = memory->create(liqtogas, nnus + 1, "kinetics/thermo:liqtogas");
  dgzero = memory->create(dgzero, atom->ntypes + 1, 2, "kinetics/thermo:dgzero");

  init_KhV();

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  } else {
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

void FixKineticsThermo::init_dG0() {
  int *ngflag = bio->ngflag;
  int *tgflag = bio->tgflag;

  for (int i = 1; i <= atom->ntypes; i++) {
    dgzero[i][0] = 0;
    dgzero[i][1] = 0;

    for (int j = 1; j <= bio->nnus; j++) {
      int flag = ngflag[j];
      dgzero[i][0] += bio->catCoeff[i][j] * bio->nuGCoeff[j][flag];
      dgzero[i][1] += bio->anabCoeff[i][j] * bio->nuGCoeff[j][flag];
    }

    int flag = tgflag[i];
    dgzero[i][1] += bio->typeGCoeff[i][flag];
  }
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_KhV() {
  int nnus = bio->nnus;
  for (int i = 1; i < nnus + 1; i++) {
    khv[i] = 0;
    liqtogas[i] = 0;
  }

  for (int i = 1; i < nnus + 1; i++) {
    char *nName;     // nutrient name
    nName = bio->nuName[i];

    if (bio->nuType[i] == 1) {
      char *lName = new char[strlen(nName)];     // corresponding liquid
      strncpy(lName, nName + 1, strlen(nName));

      for (int j = 1; j < nnus + 1; j++) {
        char *nName2;     // nutrient name
        nName2 = bio->nuName[j];
        if (strcmp(lName, nName2) == 0) {
          liqtogas[i] = j;
          if (bio->nuGCoeff[j][0] > 10000) {
            khv[j] = exp((bio->nuGCoeff[j][1] - bio->nuGCoeff[i][1]) / (-kinetics->rth * kinetics->temp));
          } else {
            khv[j] = exp((bio->nuGCoeff[j][0] - bio->nuGCoeff[i][1]) / (-kinetics->rth * kinetics->temp));
          }
          break;
        }
      }
      delete[] lName;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixKineticsThermo::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 thermodynamics
 ------------------------------------------------------------------------- */
void FixKineticsThermo::thermo() {
  int *mask = atom->mask;
  int *type = atom->type;
  int nnus = bio->nnus;

  double **DRGCat = kinetics->DRGCat;
  double **DRGAn = kinetics->DRGAn;
  double **nuR = kinetics->nuR;
  double **nuS = kinetics->nuS;
  double **gYield = kinetics->gYield;
  double ***activity = kinetics->activity;
  double **qGas = kinetics->qGas;

  init_dG0();

  double vRgT = kinetics->gvol * 1000 / (kinetics->rg * kinetics->temp);
  double vg = vol * 1000;
  double rGas, rLiq;

  if (rflag == 0) {
    for (int j = 1; j <= nnus; j++) {
      if (bio->nuType[j] == 1) {
        kinetics->nuConv[j] = true;
      }
    }
  }

  for (int i = 0; i < kinetics->bgrids; i++) {
    // gas transfer
    if (rflag == 1) {
      for (int j = 1; j <= nnus; j++) {
        nuR[j][i] = 0;
        qGas[j][i] = 0;
        if (bio->nuType[j] == 1) {
          //get corresponding liquid ID
          int liqID = liqtogas[j];
          double gasT = 0;
          if (liqID != 0) {
            if (bio->nuGCoeff[liqID][0] > 10000) {
              gasT = bio->kLa[liqID] * (activity[liqID][1][i] / khv[liqID] - activity[j][0][i]);
            } else {
              gasT = bio->kLa[liqID] * (activity[liqID][0][i] / khv[liqID] - activity[j][0][i]);
            }
            rGas = gasT;
            rLiq = -gasT * vRgT;
            // update nutrient consumption
            nuR[j][i] += rGas;
            nuR[liqID][i] += rLiq;
            //gas correction
            nuR[j][i] = nuR[j][i] * (nuS[j][i] * kinetics->rg * kinetics->temp / pressure + vg);
            //update gas concentration
            nuS[j][i] += nuR[j][i];
            qGas[j][i] = nuR[j][i] * kinetics->rg * kinetics->temp / pressure;
            if (nuS[j][i] < 0)
              nuS[j][i] = 0;
          }
        }
      }
    }

    for (int j = 1; j <= atom->ntypes; j++) {
      double rthT = kinetics->temp * kinetics->rth;

      //Gibbs free energy of the reaction
      DRGCat[j][i] = dgzero[j][0];  //catabolic energy values
      DRGAn[j][i] = dgzero[j][1] + rthT;  //anabolic energy values

      for (int k = 1; k <= nnus; k++) {
        if (bio->nuGCoeff[k][1] < 1e4) {
          double value = 0;
          int flag = bio->ngflag[k];
          if (activity[k][flag][i] == 0)
            value = 1e-20;
          else
            value = activity[k][flag][i];

          double dgr = rthT * log(value);

          DRGCat[j][i] += bio->catCoeff[j][k] * dgr;
          DRGAn[j][i] += bio->anabCoeff[j][k] * dgr;
        } else {
          error->all(FLERR, "nuGCoeff[1] is inf value");
        }
      }

      //use catabolic and anabolic energy values to derive catabolic reaction equation
      if (yflag == 1) {
        if (DRGCat[j][i] < 0) {
          gYield[j][i] = -(DRGAn[j][i] + bio->dissipation[j]) / DRGCat[j][i] + bio->eD[j];
          if (gYield[j][i] != 0)
            gYield[j][i] = 1 / gYield[j][i];
        } else {
          gYield[j][i] = 0;
        }
      }
    }
  }
}
