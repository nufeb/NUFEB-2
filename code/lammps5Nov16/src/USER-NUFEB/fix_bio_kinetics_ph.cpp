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
    if (strcmp(arg[iarg],"buffer_flag") == 0) {
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

  keq = memory->create(keq, nnus + 1, 4, "kinetics/ph:kEq");

  init_keq();
  compute_activity();
}

/* ---------------------------------------------------------------------- */

int FixKineticsPH::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

inline double sum_activity(double ***activity, double **keq, double **nus, int **nucharge, double denm, double *gsh, int w, int n, int c) {
  // not hydrated form acitivity
  activity[n][0][c] = keq[n][0] / w * nus[n][c] * gsh[2] / denm;
  // fully protonated form activity
  activity[n][1][c] = nus[n][c] * gsh[2] / denm;
  // 1st deprotonated form activity
  activity[n][2][c] = nus[n][c] * gsh[1] * keq[n][1] / denm;
  // 2nd deprotonated form activity
  activity[n][3][c] = nus[n][c] * gsh[0] * keq[n][1] * keq[n][2] / denm;
  // 3rd deprotonated form activity
  activity[n][4][c] = nus[n][c] * keq[n][1] * keq[n][2] * keq[n][3] / denm;

  double tmp[5];
  tmp[0] = nucharge[n][0] * activity[n][0][c];
  tmp[1] = nucharge[n][1] * activity[n][1][c];
  tmp[2] = nucharge[n][2] * activity[n][2][c];
  tmp[3] = nucharge[n][3] * activity[n][3][c];
  tmp[4] = nucharge[n][4] * activity[n][4][c];
  return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
}

inline void set_gsh(double *gsh, double value) {
  gsh[0] = value;
  gsh[1] = gsh[0] * gsh[0];
  gsh[2] = gsh[1] * gsh[0];
}

/* ----------------------------------------------------------------------
 buffer ph if the value is not in defined range
 ------------------------------------------------------------------------- */
void  FixKineticsPH::buffer_ph() {
  int index;
  int bgrids = kinetics->bgrids;

}


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

/* ---------------------------------------------------------------------- */

void FixKineticsPH::compute_activity() {
  int nnus = bio->nnu;
  double *sh = kinetics->sh;
  double **nus = kinetics->nus;
  double ***activity = kinetics->activity;

  double *denm = memory->create(denm, nnus + 1, "kinetics:denm");
  double gSh = pow(10, -iph);
  double gSh2 = gSh * gSh;
  double gSh3 = gSh * gSh2;

  for (int k = 1; k < nnus + 1; k++) {
    double iniNuS = bio->ini_nus[k][0];
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
    for (int j = 0; j < kinetics->bgrids; j++) {
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

 ------------------------------------------------------------------------- */

void FixKineticsPH::solve_ph() {
  if (!phflag) compute_activity();
  else dynamic_ph(0, kinetics->bgrids);
}

/* ----------------------------------------------------------------------
 compute ph field
 ------------------------------------------------------------------------- */
void FixKineticsPH::dynamic_ph(int first, int last) {
  int w = 1;

  double tol = 5e-15;
  int max_iter = 50;

  int nnus = bio->nnu;

  double *fa = memory->create(fa, kinetics->bgrids, "kinetics/ph:fa");
  double *fb = memory->create(fb, kinetics->bgrids, "kinetics/ph:fb");
  double *fun = memory->create(fun, kinetics->bgrids, "kinetics/ph:fun");
  double *df = memory->create(df, kinetics->bgrids, "kinetics/ph:df");

  double **nus = kinetics->nus;
  double temp = kinetics->temp;
  double rth = kinetics->rth;
  double ***activity = kinetics->activity;
  int **nucharge = bio->nucharge;
  double *sh = kinetics->sh;
  
  double a = 1e-14;
  double b = 1;

  for (int i = 0; i < kinetics->bgrids; i++) {
    fa[i] = a;
    fb[i] = b;
  }

  double gsh[3];
  set_gsh(gsh, a);
  for (int k = 1; k < nnus + 1; k++) {
    double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
      + keq[k][3] * keq[k][2] * keq[k][1];
    if (denm <= 0) {
      lmp->error->all(FLERR, "denm returns a zero value");
    }
#pragma ivdep
    for (int i = first; i < last; i++) {
      fa[i] += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i); 
    }
  }

  set_gsh(gsh, b);
  for (int k = 1; k < nnus + 1; k++) {
    double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
      + keq[k][3] * keq[k][2] * keq[k][1];
    if (denm <= 0) {
      lmp->error->all(FLERR, "denm returns a zero value");
    }
#pragma ivdep
    for (int i = first; i < last; i++) {
      fb[i] += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i); 
    }
  }

  bool wrong = false;
  for (int i = first; i < last; i++) {
    if (fa[i] * fb[i] > 0)
      wrong = true;
  }
  if (wrong)
    lmp->error->all(FLERR, "The sum of charges returns a wrong value");

  // Newton-Raphson method
  int ipH = 1;
  for (int i = first; i < last; i++) {
    sh[i] = pow(10, -iph);
  }

  while (ipH <= max_iter) {
    for (int i = first; i < last; i++) {
      fun[i] = 0;
      df[i] = 0;
    }

    for (int k = 1; k < nnus + 1; k++) {
#pragma ivdep
      for (int i = first; i < last; i++) {
        set_gsh(gsh, sh[i]);
        double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
          + keq[k][3] * keq[k][2] * keq[k][1];
        fun[i] += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i); 

        double ddenm = denm * denm;
        double aux = 3 * gsh[1] * (keq[k][0] / w + 1) + 2 * gsh[0] * keq[k][1] + keq[k][1] * keq[k][2];
        double tmp[5];
        tmp[0] = nucharge[k][0] * ((3 * gsh[1] * keq[k][0] * nus[k][i]) / (w * denm) - (keq[k][0] * nus[k][i] * gsh[2] * aux) / (w * ddenm));
        tmp[1] = nucharge[k][1] * ((3 * gsh[1] * nus[k][i]) / denm - (nus[k][i] * gsh[2] * aux) / ddenm);
        tmp[2] = nucharge[k][2] * ((2 * gsh[0] * keq[k][1] * nus[k][i]) / denm - (keq[k][1] * nus[k][i] * gsh[1] * aux) / ddenm);
        tmp[3] = nucharge[k][3] * ((keq[k][1] * keq[k][2] * nus[k][i]) / denm - (keq[k][1] * keq[k][2] * nus[k][i] * gsh[0] * aux) / ddenm);
        tmp[4] = nucharge[k][4] * (-(keq[k][1] * keq[k][2] * keq[k][3] * nus[k][i] * aux) / ddenm);
        df[i] += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
      }
    }

    // evaluation of the charge balance for the current Sh value, dF(Sh)
    bool flag = true;
    for (int i = first; i < last; i++) {
      double err = (fun[i] + sh[i]) / (1 + df[i]);
      sh[i] -= err;

      flag &= ((fabs(err) < 1e-14) && (fabs(fun[i] + sh[i]) < tol));
      // Checking if a valid pH was obtained
      if (flag) {
        if ((sh[i] <= 1e-14) || (sh[i] >= 1)) {
          // Fallback to false position method (regula falsi)
          ipH = 1;
          int n1 = 0;
          int n2 = 0;
          double fa_ = fa[i];
          double fb_ = fb[i];
          while (ipH < max_iter) {
            sh[i] = (fb_ * a - fa_ * b) / (fb_ - fa_);
            set_gsh(gsh, sh[i]);
            double fc_ = 0;
            for (int k = 1; k < nnus + 1; k++) {
              double denm = (1 + keq[k][0] / w) * gsh[2] + keq[k][1] * gsh[1] + keq[k][2] * keq[k][1] * gsh[0]
                + keq[k][3] * keq[k][2] * keq[k][1];
              fc_ += sum_activity(activity, keq, nus, nucharge, denm, gsh, w, k, i); 
            }
            fc_ += gsh[0];

            if (fa_ * fc_ > 0) {
              n1 += 1;
              if (n1 == 2) {
                fb_ = (fc_ / (fc_ + fa_)) * fb_;
                n1 = 0;
              }
              a = gsh[0];
              fa_ = fc_;
            } else if (fb_ * fc_ > 0) {
              n2 += 1;
              if (n2 == 2) {
                fa_ = (fc_ / (fc_ + fb_)) * fa_;
                n2 = 0;
              }
              b = gsh[0];
              fb_ = fc_;
            }

            double err1 = fabs(fc_);
            double err2 = fabs(gsh[0] - (fb_ * a - fa_ * b) / (fb_ - fa_));

            if ((err1 < tol) && (err2 < 1e-14)) {
              ipH = max_iter;
            }
            ipH++;
          }
        }
      }
    }
    if (flag) break;
    ipH++;
  }
  
  for (int k = 1; k < nnus + 1; k++) {
    if (strcmp(bio->nuname[k], "h") == 0) {
      for (int i = first; i < last; i++) {
        activity[k][1][i] = sh[i];
      }
      break;
    }
  }

  memory->destroy(fa);
  memory->destroy(fb);
  memory->destroy(fun);
  memory->destroy(df);
}

