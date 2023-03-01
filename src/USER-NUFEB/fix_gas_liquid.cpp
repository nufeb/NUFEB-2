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

#include "fix_gas_liquid.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "error.h"

#include "grid.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGasLiquid::FixGasLiquid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR,"Illegal fix nufeb/gas_liquid command");

  if (!grid->chemostat_flag)
    error->all(FLERR,"Fix reactor requires nufeb/chemostat grid style");

  iliquid = -1;
  igas = -1;
  kga = 0.0;
  h = 1.0;
  temp = 1.0;
  rg = 1.0;

  iliquid = grid->find(arg[3]);
  if (iliquid < 0)
    error->all(FLERR, "Can't find substrate(liquid) name");

  igas = grid->find(arg[4]);
  if (igas < 0)
    error->all(FLERR, "Can't find substrate(gas) name");

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "kga") == 0) {
      kga = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "h") == 0) {
      h = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      if (h <= 0)
    	error->all(FLERR, "Henry's law solubility constant (H) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "temp") == 0) {
      temp = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      if (temp <= 0)
    	error->all(FLERR, "Temperature (temp) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rg") == 0) {
      rg = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      if (rg <= 0)
    	error->all(FLERR, "Ideal gas constant must be positive");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/gas_liquid command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGasLiquid::init()
{
  if (grid->mw = nullptr)
    error->all(FLERR, "fix nufeb/gas_liquid requires molecular weight (mw) defined");
}

/* ---------------------------------------------------------------------- */

int FixGasLiquid::setmask()
{
  int mask = 0;
  mask |= CHEMISTRY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGasLiquid::chemistry_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixGasLiquid::compute()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double p_g2l, n_l2g;
  double vol = grid->cell_size * grid->cell_size * grid->cell_size;

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      p_g2l = kga * (conc[iliquid][i]/(h * grid->mw[iliquid]) - grid->bulk[igas]);
      n_l2g = -p_g2l / (rg * temp);
      // update reaction rates
      reac[igas][i] += p_g2l;
      reac[iliquid][i] += n_l2g * grid->mw[iliquid];
    }
  }
}


