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

#include "fix_diffusion_coeff.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "error.h"
#include "modify.h"
#include "grid.h"
#include "fix_diffusion_reaction.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{RATIO=0,DYNAMICS=1};

/* ---------------------------------------------------------------------- */

FixDiffusionCoeff::FixDiffusionCoeff(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix nufeb/diffusion_coeff command");

  isub = grid->find(arg[3]);
  coeff_flag = -1;
  fix_diffusion = nullptr;
  vol = 0.0;
  const_coeff = 0.0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ratio") == 0) {
      ratio = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      coeff_flag = RATIO;
      if (ratio > 1 || ratio < 0)
	error->all(FLERR, "Illegal fix nufeb/diffusion_coeff command: ratio");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dynamics") == 0) {
      coeff_flag = DYNAMICS;
      iarg += 1;
    } else {
      error->all(FLERR, "Illegal fix nufeb/diffusion_coeff command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixDiffusionCoeff::init()
{
  // find corresponding fix nufeb/diffusion/reaction for isub
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style, "nufeb/diffusion_reaction")) {
      fix_diffusion = (FixDiffusionReaction *)modify->fix[i];
      if (fix_diffusion->isub == isub) break;
      else fix_diffusion = nullptr;
    }
  }

  if (fix_diffusion == nullptr)
    error->all(FLERR, "Cannot find fix nufeb/diffusion_reaction");

  const_coeff = fix_diffusion->diff_coeff;
}

/* ---------------------------------------------------------------------- */

int FixDiffusionCoeff::setmask()
{
  int mask = 0;
  mask |= POST_PHYSICS_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionCoeff::post_physics_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixDiffusionCoeff::compute()
{
  int ncells = grid->ncells;
  double **coeff = grid->diff_coeff;
  double **dens = grid->dens;

  const_coeff = fix_diffusion->diff_coeff;
  vol = grid->cell_size * grid->cell_size * grid->cell_size;

  for (int i = 0; i < ncells; i++) {
    if (coeff_flag == RATIO) {
      if (dens[0][i] > 0) {
	    coeff[isub][i] = const_coeff * ratio;
      }
    } else if (coeff_flag == DYNAMICS) {
      coeff[isub][i] = const_coeff * (1 - (0.43 * pow(dens[0][i]/vol,0.92)) /
	  (11.19 + 0.27 * pow(dens[0][i]/vol,0.99)));
    }
  }
}


