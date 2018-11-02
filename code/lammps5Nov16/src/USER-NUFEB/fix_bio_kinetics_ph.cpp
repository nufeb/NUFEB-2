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
#include "comm.h"

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
  if (narg < 4)
    error->all(FLERR, "Not enough arguments in fix kinetics/ph command");

  //set default values
  buffer_flag = 0;
  phflag = 0;
  iph = 7.0;

  if (strcmp(arg[3], "fix") == 0)
    phflag = 0;
  else if (strcmp(arg[3], "dynamic") == 0)
    phflag = 1;
  else error->all(FLERR, "Illegal ph parameter:'fix' or 'dynamic'");

  int iarg = 4;
  while (iarg < narg){
    if (strcmp(arg[iarg],"buffer") == 0) {
      buffer_flag = force->inumeric(FLERR, arg[iarg+1]);
      if (buffer_flag != 0 && buffer_flag != 1)
        error->all(FLERR, "Illegal fix kinetics/ph command: buffer_flag");
      iarg += 2;
    } else if (strcmp(arg[iarg], "ph") == 0) {
      iph = force->numeric(FLERR, arg[iarg + 1]);
      if (iph < 0.0 || iph > 14)
        error->all(FLERR, "Illegal fix kinetics/ph command: ph");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix kinetics/ph command");
  }

}

/* ---------------------------------------------------------------------- */

FixKineticsPH::~FixKineticsPH() {
  memory->destroy(keq);
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
  int nnus = bio->nnu;

  if (bio->nnu == 0)
    error->all(FLERR, "fix_kinetics requires # of Nutrients inputs");
  else if (bio->nucharge == NULL)
    error->all(FLERR, "fix_kinetics requires Nutrient Charge inputs");

  keq = memory->create(keq, nnus + 1, 4, "kinetics/ph:keq");

  init_keq();
  compute_activity(0, kinetics->ngrids, iph);
}

/* ---------------------------------------------------------------------- */

int FixKineticsPH::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::solve_ph() {
  if (!phflag) compute_activity(0, kinetics->bgrids, iph);
  else dynamic_ph(0, kinetics->bgrids);
}


//inline double sum_activity(double ***activity, double **keq, double **nus, int **nucharge, double denm, double *gsh, int w, int n, int c) {
//  // not hydrated form acitivity
//  activity[n][0][c] = keq[n][0] / w * nus[n][c] * gsh[2] / denm;
//  // fully protonated form activity
//  activity[n][1][c] = nus[n][c] * gsh[2] / denm;
//  // 1st deprotonated form activity
//  activity[n][2][c] = nus[n][c] * gsh[1] * keq[n][1] / denm;
//  // 2nd deprotonated form activity
//  activity[n][3][c] = nus[n][c] * gsh[0] * keq[n][1] * keq[n][2] / denm;
//  // 3rd deprotonated form activity
//  activity[n][4][c] = nus[n][c] * keq[n][1] * keq[n][2] * keq[n][3] / denm;
//
//  double tmp[5];
//  tmp[0] = nucharge[n][0] * activity[n][0][c];
//  tmp[1] = nucharge[n][1] * activity[n][1][c];
//  tmp[2] = nucharge[n][2] * activity[n][2][c];
//  tmp[3] = nucharge[n][3] * activity[n][3][c];
//  tmp[4] = nucharge[n][4] * activity[n][4][c];
//  return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
//}
//
//inline void set_gsh(double *gsh, double value) {
//  gsh[0] = value;
//  gsh[1] = gsh[0] * gsh[0];
//  gsh[2] = gsh[1] * gsh[0];
//}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::init_keq() {
  // water Kj/mol
  double dG0H2O = -237.18;
  int nnus = bio->nnu;
  double **nuGCoeff = bio->nugibbs_coeff;

  for (int i = 1; i < nnus + 1; i++) {
    keq[i][0] = exp((dG0H2O + nuGCoeff[i][0] - nuGCoeff[i][1]) / (-kinetics->rth * kinetics->temp));
    for (int j = 1; j < 4; j++) {
      double coeff = 0.0;

      if (nuGCoeff[i][j + 1] > 10000) {
        coeff = j * 10001;
      } else {
        coeff = 0;
      }

      keq[i][j] = exp((nuGCoeff[i][j + 1] + coeff - nuGCoeff[i][j]) / (-kinetics->rth * kinetics->temp));
    }
  }
}

/* ----------------------------------------------------------------------
 compute nutrient form concentration, only called when fix ph is applied
 ------------------------------------------------------------------------- */

void FixKineticsPH::compute_activity(int first, int last, double iph) {
  int nnus = bio->nnu;
  double *sh = kinetics->sh;
  double **nus = kinetics->nus;
  double ***activity = kinetics->activity;

  double *denm = memory->create(denm, nnus + 1, "kinetics:denm");
  double gSh = pow(10, -iph);
  double gSh2 = gSh * gSh;
  double gSh3 = gSh * gSh2;

  for (int k = 1; k < nnus + 1; k++) {
    denm[k] = (1 + keq[k][0]) * gSh3 + keq[k][1] * gSh2 + keq[k][2] * keq[k][3] * gSh
        + keq[k][3] * keq[k][2] * keq[k][1];
    if (denm[k] == 0) {
      lmp->error->all(FLERR, "denm returns a zero value");
    }
    double tmp[5];
    tmp[0] = keq[k][0] * gSh3 / denm[k];
    tmp[1] = gSh3 / denm[k];
    tmp[2] = gSh2 * keq[k][1] / denm[k];
    tmp[3] = gSh * keq[k][1] * keq[k][2] / denm[k];
    tmp[4] = keq[k][1] * keq[k][2] * keq[k][3] / denm[k];
    bool is_hydrogen = false;
    if (strcmp(bio->nuname[k], "h") == 0) {
      is_hydrogen = true;
    }

#pragma ivdep
#pragma vector aligned
    for (int j = first; j < last; j++) {
      sh[j] = gSh;
      // not hydrated form acitivity
      activity[k][0][j] = nus[k][j] * tmp[0];
      // fully protonated form activity
      if (is_hydrogen) {
        activity[k][1][j] = gSh;
      } else {
        activity[k][1][j] = nus[k][j] * tmp[1];
      }
      // 1st deprotonated form activity
      activity[k][2][j] = nus[k][j] * tmp[2];
      // 2nd deprotonated form activity
      activity[k][3][j] = nus[k][j] * tmp[3];
      // 3rd deprotonated form activity
      activity[k][4][j] = nus[k][j] * tmp[4];
      // if(k==1)printf("act = %e, s= %e, flag = %i \n", activity[k][1][j], nus[k][j], bio->ngflag[k]);
    }
  }
  memory->destroy(denm);
}

/* ----------------------------------------------------------------------
 buffer ph if the value is not in defined range
 ------------------------------------------------------------------------- */
void  FixKineticsPH::buffer_ph() {
  int nnus = bio->nnu;

  if (bio->find_nuid("na") < 0 || !bio->find_nuid("cl") < 0)
    error->all(FLERR, "buffer ph requires nutreint 'na' and 'cl'");

  int grid;
  double ph_unbuffer;
  int **nucharge = bio->nucharge;
  double ***activity = kinetics->activity;

  // always take the last grid
  grid = kinetics->ngrids - 1;
  dynamic_ph(grid, grid+1);
  ph_unbuffer = -log10(kinetics->sh[grid]);

  if (ph_unbuffer < 6.5 || ph_unbuffer > 9) {
    double minus = 0;
    double plus = 0;
    compute_activity(grid, grid+1, 7);

    for (int nu = 1; nu <= nnus ; nu++){
      for (int i = 0; i < 5; i++) {
        double diff, act, chr;
        act = activity[nu][i][grid];
        chr = nucharge[nu][i];

        if (act > 10000 || chr > 10000) continue;

        diff = act * chr;

        if (diff > 0) plus += diff;
        else if (diff < 0) minus -= diff;
      }
    }

    kinetics->nubs[bio->find_nuid("na")] += minus;
    kinetics->nubs[bio->find_nuid("cl")] += plus + kinetics->sh[grid];
  }
}

/* ----------------------------------------------------------------------
 compute ph field
 ------------------------------------------------------------------------- */
void FixKineticsPH::dynamic_ph(int first, int last) {
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
  int **nucharge = bio->nucharge;

  for (int i = first; i < last; i++) {
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
    gsh = pow(10, -iph);

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

