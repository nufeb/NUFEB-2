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

#include "fix_kinetics_ph.h"

#include <math.h>
#include <string.h>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <cfloat>

#include "atom.h"
#include "bio.h"
#include "domain.h"
#include "error.h"
#include "fix_kinetics.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "memory.h"
#include "modify.h"
#include "pointers.h"
#include "variable.h"
#include "update.h"


using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsPH::FixKineticsPH(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics/ph command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/ph command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsPH::~FixKineticsPH()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  memory->destroy(kEq);
  memory->destroy(Sh);
}

/* ---------------------------------------------------------------------- */

int FixKineticsPH::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::init()
{
  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics/ph does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics/ph is invalid style");
  }

  ph = input->variable->compute_equal(ivar[0]);

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for kinetics/monod styles");

  ntypes = atom->ntypes;

  bio = kinetics->bio;
  nnus = bio->nnus;
  nuGCoeff = bio->nuGCoeff;
  typeGCoeff = bio->typeGCoeff;
  nuChr = bio->nuChr;

  temp = kinetics->temp;
  rth = kinetics->rth;
  activity = kinetics->activity;

  Sh = memory->create(Sh,kinetics->ngrids,"kinetics/ph:Sh");
  kEq = memory->create(kEq,nnus+1,4,"kinetics/ph:nuKeq");
  init_keq();
}

/* ---------------------------------------------------------------------- */

void FixKineticsPH::init_keq()
{
  // water Kj/mol
  double dG0H2O = -237.18;

  for (int i = 1; i < nnus+1; i++) {
    kEq[i][0] = exp((dG0H2O + nuGCoeff[i][0] - nuGCoeff[i][1]) / (-rth * temp));
    for (int j = 1; j < 4; j++) {
      double coeff = 0.0;

      if (nuGCoeff[i][j+1] > 10000) {
        coeff = j * 10001;
      } else {
        coeff = 0;
      }

      kEq[i][j] = exp((nuGCoeff[i][j+1] + coeff - nuGCoeff[i][j]) / (-rth * temp));
    }
  }
}


/* ---------------------------------------------------------------------- */

void FixKineticsPH::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  solve_ph();
}

/* ----------------------------------------------------------------------
  ph calculation
------------------------------------------------------------------------- */
void FixKineticsPH::solve_ph()
{
  int w = 1; // why divide by w???
  double a = 1e-14;
  int b = 1;

  int ipH = 1;
  double tol = 5e-15;
  int maxIter = 20;
  double fun;
  double dF;
  double fa, fb, fc;
  double err, err1, err2;
  double sumActivity;

  double *denm = memory->create(denm,nnus+1,"kinetics/ph:denm");
  double *dDenm = memory->create(dDenm,nnus+1,"kinetics/ph:denm");
  double *aux = memory->create(aux,nnus+1,"kinetics/ph:aux");

  nuS = kinetics->nuS;

  for (int i = 0; i < kinetics->ngrids; i++) {
    double gSh; // grid Sh
    for (int j = 0; j < 2; j++) {
      sumActivity = 0.0;
      if (j == 0) gSh = a;
      else gSh = b;
      for (int k = 1; k < nnus+1; k++) {
        denm[k] = (1 + kEq[k][0]/w) * gSh * gSh * gSh + kEq[k][1] * gSh * gSh + kEq[k][2] * kEq[k][3] * gSh + kEq[k][3] * kEq[k][2] * kEq[k][1];
        if (denm[k] == 0) {
          lmp->error->all(FLERR,"denm returns a zero value");
        }
        // not hydrated form acitivity
        activity[k][0] = kEq[k][0]/w * nuS[k][i] * gSh * gSh * gSh / denm[k];
        // fully protonated form activity
        activity[k][1] = nuS[k][i] * gSh * gSh * gSh / denm[k];
        // 1st deprotonated form activity
        activity[k][2] = nuS[k][i] * gSh * gSh * kEq[k][1] / denm[k];
        // 2nd deprotonated form activity
        activity[k][3] = nuS[k][i] * gSh * kEq[k][1] * kEq[k][2] / denm[k];
        // 3rd deprotonated form activity
        activity[k][4] = nuS[k][i] * kEq[k][1] * kEq[k][2] * kEq[k][3] / denm[k];

        sumActivity += nuChr[k][0] * activity[k][0];
        sumActivity += nuChr[k][1] * activity[k][1];
        sumActivity += nuChr[k][2] * activity[k][2];
        sumActivity += nuChr[k][3] * activity[k][3];
        sumActivity += nuChr[k][4] * activity[k][4];
      }

      if (j == 0) fa = gSh + sumActivity;
      else fb = gSh + sumActivity;
    }

    double ff = fa * fb;
    if (ff > 0) {
      lmp->error->all(FLERR,"The sum of charges returns a wrong value");
    }

    //Newton-Raphson method
    gSh = pow(10, -ph);

    while (ipH <= maxIter) {
      sumActivity = 0.0;
      for (int k = 1; k < nnus+1; k++) {
        denm[k] = (1 + kEq[k][0]/w) * gSh * gSh * gSh + kEq[k][1] * gSh * gSh + kEq[k][2] * kEq[k][1] * gSh + kEq[k][3] * kEq[k][2] * kEq[k][1] ;

        activity[k][0] = kEq[k][0]/w * nuS[k][i] * gSh * gSh * gSh / denm[k];
        activity[k][1] = nuS[k][i] * gSh * gSh * gSh / denm[k];
        activity[k][2] = nuS[k][i] * gSh * gSh * kEq[k][1] / denm[k];
        activity[k][3] = nuS[k][i] * gSh * kEq[k][1] * kEq[k][2] / denm[k];
        activity[k][4] = nuS[k][i] * kEq[k][1] * kEq[k][2] * kEq[k][3] / denm[k] ;

        sumActivity += nuChr[k][0] * activity[k][0];
        sumActivity += nuChr[k][1] * activity[k][1];
        sumActivity += nuChr[k][2] * activity[k][2];
        sumActivity += nuChr[k][3] * activity[k][3];
        sumActivity += nuChr[k][4] * activity[k][4];
      }
      // evaluation of the charge balance for the current Sh value, F(Sh)
      fun = gSh + sumActivity;

      sumActivity = 0.0;
      for (int k = 1; k < nnus+1; k++) {
        dDenm[k] = pow(denm[k], 2);
        aux[k] = 3 * gSh * gSh * (kEq[k][0]/w + 1) + 2 * gSh * kEq[k][1] + kEq[k][1] * kEq[k][2];

        sumActivity += nuChr[k][0] * ((3 * gSh * gSh * kEq[k][0] * nuS[k][i]) / (w * denm[k]) - (kEq[k][0] * nuS[k][i] * gSh * gSh * gSh * aux[k]) / (w * dDenm[k]));
        sumActivity += nuChr[k][1] * ((3 * gSh * gSh * nuS[k][i]) / denm[k] - (nuS[k][i] * gSh * gSh * gSh * aux[k]) / dDenm[k]);
        sumActivity += nuChr[k][2] * ((2 * gSh * kEq[k][1] * nuS[k][i]) / denm[k] - (kEq[k][1] * nuS[k][i] * gSh * gSh * aux[k]) / dDenm[k]);
        sumActivity += nuChr[k][3] * ((kEq[k][1] * kEq[k][2] * nuS[k][i]) / denm[k] - (kEq[k][1] * kEq[k][2] * nuS[k][i] * gSh * aux[k]) / dDenm[k]);
        sumActivity += nuChr[k][4] * (-(kEq[k][1] * kEq[k][2] * kEq[k][3] * nuS[k][i] * aux[k]) / dDenm[k]);
      }
      // evaluation of the charge balance for the current Sh value, dF(Sh)
      dF = 1 + sumActivity;
      err = fun/dF;
      gSh = gSh - err;

      if ((abs(err) < 1e-14) && (abs(fun) < tol)) {
        // Checking if a valid pH was obtained
        if ((gSh > 1e-14) && (gSh < 1)) {
          ipH = maxIter;
        } else {
          //Counter of convergence
          ipH = 1;
          maxIter = 50;
          int n1 = 0;
          int n2 = 0;
          while (ipH < maxIter) {
            gSh = (fa * a - fa * b) / (fb - fa);
            sumActivity = 0.0;
            for (int k = 1; k < nnus+1; k++) {
              denm[k] = (1 + kEq[k][0]/w) * gSh * gSh * gSh + kEq[k][1] * gSh * gSh + kEq[k][2] * kEq[k][1] * gSh + kEq[k][3] * kEq[k][2] * kEq[k][1];

              activity[k][0] = kEq[k][0]/w * nuS[k][i] * gSh * gSh * gSh / denm[k];
              activity[k][1] = nuS[k][i] * gSh * gSh * gSh / denm[k];
              activity[k][2] = nuS[k][i] * gSh * gSh * kEq[k][1] / denm[k];
              activity[k][3] = nuS[k][i] * gSh * kEq[k][1] * kEq[k][2] / denm[k];
              activity[k][4] = nuS[k][i] * kEq[k][1] * kEq[k][2] * kEq[k][3] / denm[k];

              sumActivity += nuChr[k][0] * activity[k][0];
              sumActivity += nuChr[k][1] * activity[k][1];
              sumActivity += nuChr[k][2] * activity[k][2];
              sumActivity += nuChr[k][3] * activity[k][3];
              sumActivity += nuChr[k][4] * activity[k][4];
            }
            fc = gSh + sumActivity;
            if (fa * fc > 0) {
              n1 += 1;
              if (n1 == 2) {
                fb = (fc / (fc + fa)) * fb;
                n1 = 0;
              }
              a = gSh;
              fa = fc;
            } else if (fb * fc > 0) {
              n2 += 1;
              if (n2 == 2) {
                fa = (fc / (fc + fb)) * fa;
                n2 = 0;
              }
              b = gSh;
              fb = fc;
            }
            err1 = abs(fc);
            err2 = abs(gSh - (fb * a - fa * b) / (fb - fa));

            if ((err1 < tol) && (err2 < 1e-14)) {
              ipH = maxIter;
            }
            ipH += 1;
          }
        }
      }
      ipH += ipH;
    }
    Sh[i] = gSh;
  }

  memory->destroy(dDenm);
  memory->destroy(denm);
  memory->destroy(aux);
}

/* ----------------------------------------------------------------------
  output energy to data file
------------------------------------------------------------------------- */

void FixKineticsPH::output_data(){


}
