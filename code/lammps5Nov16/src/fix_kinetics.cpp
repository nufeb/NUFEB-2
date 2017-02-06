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

#include "fix_kinetics.h"

#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "atom_vec_bio.h"
#include "bio.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "pointers.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKinetics::FixKinetics(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 10) error->all(FLERR,"Not enough arguments in fix kinetics command");

  var = new char*[4];
  ivar = new int[4];

  nx = atoi(arg[3]);
  ny = atoi(arg[4]);
  nz = atoi(arg[5]);

  for (int i = 0; i < 4; i++) {
    int n = strlen(&arg[6+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[6+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixKinetics::~FixKinetics()
{
  int i;
  for (i = 0; i < 4; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;

  memory->destroy(metCoeff);
  memory->destroy(iyield);
  memory->destroy(activity);
  memory->destroy(nuR);
  memory->destroy(nuS);
  memory->destroy(nuGas);
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
  for (int n = 0; n < 4; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics is invalid style");
  }

  temp = input->variable->compute_equal(ivar[0]);
  rth = input->variable->compute_equal(ivar[1]);
  gVol = input->variable->compute_equal(ivar[2]);
  gasTrans = input->variable->compute_equal(ivar[3]);

  bio = avec->bio;

  ngrids = nx * ny * nz;

  int nnus = bio->nnus;
  int ntypes = atom->ntypes;

  nuS = memory->create(nuS,nnus+1, ngrids, "kinetics:nuS");
  nuR = memory->create(nuR,nnus+1, ngrids, "kinetics:nuR");
  nuGas = memory->create(nuGas,nnus+1, ngrids, "kinetics:nuGas");
  metCoeff = memory->create(metCoeff,ntypes+1,nnus+1,"kinetic:metCoeff");
  iyield = memory->create(iyield,ntypes+1,ngrids,"kinetic:iyield");
  activity = memory->create(activity,nnus+1,5,"kinetics/ph:activity");

  //initialize inlet concentration, consumption
  for (int i = 1; i <= nnus; i++) {
    for (int j = 0; j < ngrids; j++) {
      nuS[i][j] = bio->iniS[i][0];
      nuR[i][j] = 0;
    }
  }
}


