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
#include "atom_vec_bio.h"
#include "bio.h"
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
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKinetics::FixKinetics(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 7) error->all(FLERR,"Not enough arguments in fix kinetics command");

  var = new char*[1];
  ivar = new int[1];

  nx = atoi(arg[3]);
  ny = atoi(arg[4]);
  nz = atoi(arg[5]);

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[6+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[6+i][2]);
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
  memory->destroy(iyield);
  memory->destroy(nuR);
  memory->destroy(nuS);
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
  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix kinetics is invalid style");
  }

  bio = avec->bio;

  ngrids = nx * ny * nz;
  temp = input->variable->compute_equal(ivar[0]);

  int nnus = bio->nnus;
  int ntypes = atom->ntypes;
  nuS = memory->create(nuS,nnus+1, ngrids, "kinetics:nuS");
  nuR = memory->create(nuR,nnus+1, ngrids, "kinetics:nuR");
  metCoeff = memory->create(metCoeff,ntypes+1,nnus+1,"kinetic:metCoeff");
  iyield = memory->create(iyield,ntypes+1,ngrids,"kinetic:iyield");
}


