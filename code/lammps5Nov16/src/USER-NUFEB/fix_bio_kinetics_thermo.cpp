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

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "yield") == 0) {
      if (strcmp(arg[iarg + 1], "fix") == 0)
        yflag = 0;
      else if (strcmp(arg[iarg + 1], "unfix") == 0)
        yflag = 1;
      else
        error->all(FLERR, "Illegal yield parameter:fix or unfix");
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor") == 0) {
      if (strcmp(arg[iarg + 1], "close") == 0)
        rflag = 1;
      else if (strcmp(arg[iarg + 1], "open") == 0)
        rflag = 0;
      else
        error->all(FLERR, "Illegal reactor parameter:open or close");
      iarg += 2;
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

  if (bio->cata_coeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Catabolism Coeffs input");
  else if (bio->anab_coeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Anabolism Coeffs input");
  else if (bio->nugibbs_coeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Nutrient Energy input");
  else if (bio->tgibbs_coeff == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Type Energy input");
  else if (yflag == 1 && bio->dissipation == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires Dissipation inputs for unfix yield");
  else if (rflag == 1 && bio->kla == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires KLa input for closed reactor");
  else if (yflag == 1 && bio->edoner == NULL)
    error->all(FLERR, "fix_kinetics/thermo requires eD input for unfix yield");

  int nnus = bio->nnu;

  khv = memory->create(khv, nnus + 1, "kinetics/thermo:khv");
  liqtogas = memory->create(liqtogas, nnus + 1, "kinetics/thermo:liqtogas");
  dgzero = memory->create(dgzero, atom->ntypes + 1, 2, "kinetics/thermo:dgzero");

  init_KhV();
  init_dG0();

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
  int *tgflag = bio->tgibbs_flag;

  for (int i = 1; i <= atom->ntypes; i++) {
    dgzero[i][0] = 0;
    dgzero[i][1] = 0;

    for (int j = 1; j <= bio->nnu; j++) {
      int flag = ngflag[j];
      dgzero[i][0] += bio->cata_coeff[i][j] * bio->nugibbs_coeff[j][flag];
      dgzero[i][1] += bio->anab_coeff[i][j] * bio->nugibbs_coeff[j][flag];
    }

    int flag = tgflag[i];
    dgzero[i][1] += bio->tgibbs_coeff[i][flag];
  }
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_KhV() {
  int nnus = bio->nnu;
  for (int i = 1; i < nnus + 1; i++) {
    khv[i] = 0;
    liqtogas[i] = 0;
  }

  for (int i = 1; i < nnus + 1; i++) {
    char *nName;     // nutrient name
    nName = bio->nuname[i];

    if (bio->nustate[i] == 1) {
      char *lName = new char[strlen(nName)];     // corresponding liquid
      strncpy(lName, nName + 1, strlen(nName));

      for (int j = 1; j < nnus + 1; j++) {
        char *nName2;     // nutrient name
        nName2 = bio->nuname[j];
        if (strcmp(lName, nName2) == 0) {
          liqtogas[i] = j;

          if (bio->nugibbs_coeff[j][0] > 10000) {
            khv[j] = exp((bio->nugibbs_coeff[j][1] - bio->nugibbs_coeff[i][1]) / (-kinetics->rth * kinetics->temp));
          } else {
            khv[j] = exp((bio->nugibbs_coeff[j][0] - bio->nugibbs_coeff[i][1]) / (-kinetics->rth * kinetics->temp));
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
void FixKineticsThermo::thermo(double dt) {
  int *mask = atom->mask;
  int *type = atom->type;
  int nnus = bio->nnu;

  double **gibbs_cata = kinetics->gibbs_cata;
  double **gibbs_anab = kinetics->gibbs_anab;
  double **nur = kinetics->nur;
  double **nus = kinetics->nus;
  double **grid_yield = kinetics->grid_yield;
  double ***activity = kinetics->activity;

  double vRgT = kinetics->gvol * 1000 / (kinetics->rg * kinetics->temp);
  double vg = vol * 1000;
  double rGas, rLiq;

  if (rflag == 0) {
    for (int j = 1; j <= nnus; j++) {
      if (bio->nustate[j] == 1)
        kinetics->nuconv[j] = true;
    }
  }

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    // gas liquid transfer
    if (rflag == 1) {
      for (int nu = 1; nu <= nnus; nu++) {
        if (bio->nustate[nu] != 1) continue;

        //get corresponding liquid ID
        int liqID = liqtogas[nu];
        double gasT = 0;

        if (!liqID) continue;
        if (bio->nugibbs_coeff[liqID][0] > 10000) {
          gasT = bio->kla[liqID] * (activity[nu][0][grid] - activity[liqID][1][grid] / khv[liqID]);
        } else {
          gasT = bio->kla[liqID] * (activity[nu][0][grid] - activity[liqID][0][grid] / khv[liqID]);
        }

        rGas = -gasT;
        rLiq = gasT * vRgT;
        // update nutrient consumption
        nur[nu][grid] += rGas;
        nur[liqID][grid] += rLiq;

        //update gas concentration
        nus[nu][grid] += nur[nu][grid] * dt;

        if (nus[nu][grid] < 0) {
          nus[nu][grid] = 1e-20;
        }
      }
    }

    // calculate metablic energy
    for (int i = 1; i <= atom->ntypes; i++) {
      double rthT = kinetics->temp * kinetics->rth;

      //Gibbs free energy of the reaction
      gibbs_cata[i][grid] = dgzero[i][0];  //catabolic energy values
      gibbs_anab[i][grid] = dgzero[i][1] + rthT;  //anabolic energy values

      for (int nu = 1; nu <= nnus; nu++) {
        if (bio->nugibbs_coeff[nu][1] >= 1e4)
          error->all(FLERR, "nuGCoeff[1] is inf value");

        double act = 0;
        int flag = bio->ngflag[nu];
        if (activity[nu][flag][grid] == 0)
          act = 1e-20;
        else
          act = activity[nu][flag][grid];

        double dgr = rthT * log(act);

        gibbs_cata[i][grid] += bio->cata_coeff[i][nu] * dgr;
        gibbs_anab[i][grid] += bio->anab_coeff[i][nu] * dgr;
      }

      //use catabolic and anabolic energy values to derive catabolic reaction equation
      if (!yflag) continue;
      if (gibbs_cata[i][grid] < 0) {
        grid_yield[i][grid] = -(gibbs_anab[i][grid] + bio->dissipation[i]) / gibbs_cata[i][grid] + bio->edoner[i];
        if (grid_yield[i][grid] != 0)
          grid_yield[i][grid] = 1 / grid_yield[i][grid];
      } else {
        grid_yield[i][grid] = 0;
      }

    }
  }
}

