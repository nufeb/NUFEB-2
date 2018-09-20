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

#include "fix_bio_kinetics_ph.h"

#include <math.h>
#include <string.h>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <cfloat>

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

FixKineticsPH::FixKineticsPH(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  if (narg != 3)
    error->all(FLERR, "Not enough arguments in fix kinetics/ph command");
}

/* ---------------------------------------------------------------------- */

FixKineticsPH::~FixKineticsPH() {
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::init() {
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
    lmp->error->all(FLERR, "The fix kinetics command is required for kinetics/ph styles");

  bio = kinetics->bio;

  if (bio->nnu == 0)
    error->all(FLERR, "fix_kinetics requires # of Nutrients inputs");
  else if (bio->nucharge == NULL)
    error->all(FLERR, "fix_kinetics requires Nutrient Charge inputs");
}

/* ---------------------------------------------------------------------- */

int FixKineticsPH::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 ph calculation
 ------------------------------------------------------------------------- */
void FixKineticsPH::solve_ph() {
  int w = 1;

  double tol = 5e-15;
  int max_iter = 50;

  int nnus = bio->nnu;

  double *denm = memory->create(denm, nnus + 1, "kinetics/ph:denm");
  double *ddenm = memory->create(ddenm, nnus + 1, "kinetics/ph:denm");
  double *aux = memory->create(aux, nnus + 1, "kinetics/ph:aux");

  double **nus = kinetics->nus;
  double temp = kinetics->temp;
  double rth = kinetics->rth;
  double ***activity = kinetics->activity;
  double **keq = kinetics->keq;
  int **nucharge = bio->nucharge;

  for (int i = 0; i < kinetics->bgrids; i++) {
    double a = 1e-14;
    double b = 1;
    int ipH = 1;

    double fun;
    double df;
    double fa, fb, fc;
    double err, err1, err2;
    double sum_activity;
    double gsh; // grid Sh

    for (int j = 0; j < 2; j++) {
      sum_activity = 0.0;
      if (j == 0)
        gsh = a;
      else
        gsh = b;
      for (int k = 1; k < nnus + 1; k++) {
        denm[k] = (1 + keq[k][0] / w) * gsh * gsh * gsh + keq[k][1] * gsh * gsh + keq[k][2] * keq[k][1] * gsh
            + keq[k][3] * keq[k][2] * keq[k][1];
        if (denm[k] <= 0) {
          lmp->error->all(FLERR, "denm returns a zero value");
        }
        // not hydrated form acitivity
        activity[k][0][i] = keq[k][0] / w * nus[k][i] * gsh * gsh * gsh / denm[k];
        // fully protonated form activity
        activity[k][1][i] = nus[k][i] * gsh * gsh * gsh / denm[k];
        // 1st deprotonated form activity
        activity[k][2][i] = nus[k][i] * gsh * gsh * keq[k][1] / denm[k];
        // 2nd deprotonated form activity
        activity[k][3][i] = nus[k][i] * gsh * keq[k][1] * keq[k][2] / denm[k];
        // 3rd deprotonated form activity
        activity[k][4][i] = nus[k][i] * keq[k][1] * keq[k][2] * keq[k][3] / denm[k];

        sum_activity += nucharge[k][0] * activity[k][0][i];
        sum_activity += nucharge[k][1] * activity[k][1][i];
        sum_activity += nucharge[k][2] * activity[k][2][i];
        sum_activity += nucharge[k][3] * activity[k][3][i];
        sum_activity += nucharge[k][4] * activity[k][4][i];
      }

      if (j == 0)
        fa = gsh + sum_activity;
      else
        fb = gsh + sum_activity;
    }

    double ff = fa * fb;

    if (ff > 0) {
      lmp->error->all(FLERR, "The sum of charges returns a wrong value");
    }

    //Newton-Raphson method
    gsh = pow(10, -kinetics->iph);

    while (ipH <= max_iter) {
      sum_activity = 0.0;
      for (int k = 1; k < nnus + 1; k++) {
        denm[k] = (1 + keq[k][0] / w) * gsh * gsh * gsh + keq[k][1] * gsh * gsh + keq[k][2] * keq[k][1] * gsh
            + keq[k][3] * keq[k][2] * keq[k][1];

        activity[k][0][i] = keq[k][0] / w * nus[k][i] * gsh * gsh * gsh / denm[k];
        activity[k][1][i] = nus[k][i] * gsh * gsh * gsh / denm[k];
        activity[k][2][i] = nus[k][i] * gsh * gsh * keq[k][1] / denm[k];
        activity[k][3][i] = nus[k][i] * gsh * keq[k][1] * keq[k][2] / denm[k];
        activity[k][4][i] = nus[k][i] * keq[k][1] * keq[k][2] * keq[k][3] / denm[k];

        sum_activity += nucharge[k][0] * activity[k][0][i];
        sum_activity += nucharge[k][1] * activity[k][1][i];
        sum_activity += nucharge[k][2] * activity[k][2][i];
        sum_activity += nucharge[k][3] * activity[k][3][i];
        sum_activity += nucharge[k][4] * activity[k][4][i];
      }
      // evaluation of the charge balance for the current Sh value, F(Sh)
      fun = gsh + sum_activity;

      sum_activity = 0.0;
      for (int k = 1; k < nnus + 1; k++) {
        ddenm[k] = pow(denm[k], 2);
        aux[k] = 3 * gsh * gsh * (keq[k][0] / w + 1) + 2 * gsh * keq[k][1] + keq[k][1] * keq[k][2];

        sum_activity += nucharge[k][0] * ((3 * gsh * gsh * keq[k][0] * nus[k][i]) / (w * denm[k]) - (keq[k][0] * nus[k][i] * gsh * gsh * gsh * aux[k]) / (w * ddenm[k]));
        sum_activity += nucharge[k][1] * ((3 * gsh * gsh * nus[k][i]) / denm[k] - (nus[k][i] * gsh * gsh * gsh * aux[k]) / ddenm[k]);
        sum_activity += nucharge[k][2] * ((2 * gsh * keq[k][1] * nus[k][i]) / denm[k] - (keq[k][1] * nus[k][i] * gsh * gsh * aux[k]) / ddenm[k]);
        sum_activity += nucharge[k][3] * ((keq[k][1] * keq[k][2] * nus[k][i]) / denm[k] - (keq[k][1] * keq[k][2] * nus[k][i] * gsh * aux[k]) / ddenm[k]);
        sum_activity += nucharge[k][4] * (-(keq[k][1] * keq[k][2] * keq[k][3] * nus[k][i] * aux[k]) / ddenm[k]);
      }
      // evaluation of the charge balance for the current Sh value, dF(Sh)
      df = 1 + sum_activity;
      err = fun / df;
      gsh = gsh - err;

      if ((fabs(err) < 1e-14) && (fabs(fun) < tol)) {
        // Checking if a valid pH was obtained
        if ((gsh > 1e-14) && (gsh < 1)) {
          ipH = max_iter;
        } else {
          //Counter of convergence
          ipH = 1;
          max_iter = 50;
          int n1 = 0;
          int n2 = 0;
          while (ipH < max_iter) {
            gsh = (fb * a - fa * b) / (fb - fa);
            sum_activity = 0.0;

            for (int k = 1; k < nnus + 1; k++) {
              denm[k] = (1 + keq[k][0] / w) * gsh * gsh * gsh + keq[k][1] * gsh * gsh + keq[k][2] * keq[k][1] * gsh
                  + keq[k][3] * keq[k][2] * keq[k][1];

              activity[k][0][i] = keq[k][0] / w * nus[k][i] * gsh * gsh * gsh / denm[k];
              activity[k][1][i] = nus[k][i] * gsh * gsh * gsh / denm[k];
              activity[k][2][i] = nus[k][i] * gsh * gsh * keq[k][1] / denm[k];
              activity[k][3][i] = nus[k][i] * gsh * keq[k][1] * keq[k][2] / denm[k];
              activity[k][4][i] = nus[k][i] * keq[k][1] * keq[k][2] * keq[k][3] / denm[k];

              sum_activity += nucharge[k][0] * activity[k][0][i];
              sum_activity += nucharge[k][1] * activity[k][1][i];
              sum_activity += nucharge[k][2] * activity[k][2][i];
              sum_activity += nucharge[k][3] * activity[k][3][i];
              sum_activity += nucharge[k][4] * activity[k][4][i];
            }
            fc = gsh + sum_activity;

            if (fa * fc > 0) {
              n1 += 1;
              if (n1 == 2) {
                fb = (fc / (fc + fa)) * fb;
                n1 = 0;
              }
              a = gsh;
              fa = fc;
            } else if (fb * fc > 0) {
              n2 += 1;
              if (n2 == 2) {
                fa = (fc / (fc + fb)) * fa;
                n2 = 0;
              }
              b = gsh;
              fb = fc;
            }

            err1 = fabs(fc);
            err2 = fabs(gsh - (fb * a - fa * b) / (fb - fa));

            if ((err1 < tol) && (err2 < 1e-14)) {
              ipH = max_iter;
            }
            ipH++;
          }
        }
      }
      ipH++;
    }

    kinetics->sh[i] = gsh;

    for (int k = 1; k < nnus + 1; k++) {
      if (strcmp(bio->nuname[k], "h") == 0) {
        activity[k][1][i] = gsh;
        break;
      }
    }
  }

  memory->destroy(ddenm);
  memory->destroy(denm);
  memory->destroy(aux);
}

