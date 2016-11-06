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

#include <fix_kinetics.h>
#include <fix_diffusion.h>
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
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

FixKinetics::FixKinetics(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Not enough arguments in fix kinetics command");

  nx = ny = nz = 1;
  nnus = 0;
  ntypes = 0;
  temp = 0.0;

  catCoeff = NULL;
  anabCoeff = NULL;
  metCoeff = NULL;
  yield = NULL;

  nuS = NULL;
  nuR = NULL;
  nuG = NULL;

  diffusion = NULL;

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixKinetics::~FixKinetics()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  memory->destroy(metCoeff);
}

/* ---------------------------------------------------------------------- */

int FixKinetics::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init()
{
  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"diffusion") == 0) {
      diffusion = static_cast<FixDiffusion *>(lmp->modify->fix[j]);
      break;
    }
  }

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics is invalid style");
  }

  temp = input->variable->compute_equal(ivar[0]);

  if (diffusion != NULL) {
    nx = diffusion->nx;
    ny = diffusion->ny;
    nz = diffusion->nz;
  }

  ngrids = nx * ny * nz;
  nnus = atom->nNutrients;
  ntypes = atom->ntypes;
  catCoeff = atom->catCoeff;
  anabCoeff = atom->anabCoeff;
  nuGCoeff = atom->nuGCoeff;
  typeGCoeff = atom->typeGCoeff;
  yield = atom->yield;

  nuS = atom->nuS;
  nuR = atom->nuR;
  nuG = atom->nuGCoeff;

  metCoeff = memory->create(metCoeff,ntypes+1,nnus+1,"kinetic:metCoeff");

}


