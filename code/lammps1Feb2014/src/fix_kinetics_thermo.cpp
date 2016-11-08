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

#include "fix_kinetics_thermo.h"
#include <fix_kinetics.h>
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

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsThermo::FixKineticsThermo(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics/thermo command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/thermo command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsThermo::~FixKineticsThermo()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  memory->destroy(dG0);
  memory->destroy(iyield);
//  memory->destroy(dGrxn);
}

/* ---------------------------------------------------------------------- */

int FixKineticsThermo::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::init()
{

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics is invalid style");
  }

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

  rth = input->variable->compute_equal(ivar[0]);

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  ngrids = nx * ny * nz;
  nnus = kinetics->nnus;
  ntypes = kinetics->ntypes;
  catCoeff = kinetics->catCoeff;
  anabCoeff = kinetics->anabCoeff;
  metCoeff = kinetics->metCoeff;
  nuGCoeff = kinetics->nuGCoeff;
  typeGCoeff = kinetics->typeGCoeff;
  yield = kinetics->yield;
  temp = kinetics->temp;
  diss = kinetics->dissipation;

  iyield = memory->create(iyield,ntypes+1,ngrids,"kinetic/thermo:iyield");
  //dGrxn = memory->create(dGrxn,ntypes+1,2,ngrids,"kinetics/thermo:dGrxn");
  init_dG0();
}

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

void FixKineticsThermo::init_dG0()
{
  dG0 = memory->create(dG0,ntypes+1,2,"kinetics/thermo:dG0");

  for (int i = 1; i <= ntypes; i++) {
    dG0[i][0] = 0;
    dG0[i][1] = 0;

    for (int j = 1; j <= nnus; j++) {
      if (nuGCoeff[j][1] < 1e4) {
        dG0[i][0] += catCoeff[i][j] * nuGCoeff[j][1];
        dG0[i][1] += anabCoeff[i][j] * nuGCoeff[j][1];
      }
    }

    if (typeGCoeff[i][1] < 1e4) {
      dG0[i][1] += typeGCoeff[i][1];
    }
//
//    printf("dG0[%i][0] = %e \n", i, dG0[i][0]);
//    printf("dG0[%i][1] = %e \n", i, dG0[i][1]);
  }
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  thermo();
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

  nuS = kinetics->nuS;

  for (int i = 0; i < ngrids; i++) {
    for (int j = 1; j <= ntypes; j++) {
      //Gibbs free energy of the reaction
      double catG = dG0[j][0];  //catabolic energy values
      double anaG = dG0[j][1];  //anabolic energy values

      for (int k = 1; k <= nnus; k++) {
        double x = temp * rth * log(nuS[k][i]);

        if (nuGCoeff[k][1] < 1e4) {
          catG += catCoeff[j][k] * x;
          anaG += anabCoeff[j][k] * x;
        }
      }

//      printf("dGrxn[%i][%i][0] = %e \n", i,j, catG);
//      printf("dGrxn[%i][%i][1] = %e \n", i,j, anaG);

      //use catabolic and anabolic energy values to derive catabolic reaction equation
      iyield[j][i] = - anaG * (- (anaG + diss[j]) / catG);
//      printf("iyield = %e \n",  iyield[j][i]);
    }
  }
}

