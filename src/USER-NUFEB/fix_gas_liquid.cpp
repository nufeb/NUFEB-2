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

#include <cstdio>
#include <cstring>
#include <cmath>
#include "error.h"

#include "fix_reactor_gas_liquid.h"
#include "grid.h"
#include "grid_masks.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReactorGasLiquid::FixReactorGasLiquid(LAMMPS *lmp, int narg, char **arg) :
  FixReactor(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR,"Illegal fix nufeb/gas_liquid command");

  iliquid = -1;
  igas = -1;

  kga = 0.0;
  h = 1.0;
  temp = 1.0;
  reactor_vhead = 1.0;
  reactor_pres = 1.0;
  mw = 1.0;
  rg = 1.0;

  rtotal = 0.0;

  iliquid = grid->find(arg[3]);
  if (iliquid < 0)
    error->all(FLERR, "Can't find substrate(liquid) name");

  igas = grid->find(arg[4]);
  if (igas < 0)
    error->all(FLERR, "Can't find substrate(gas) name");

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "kla") == 0) {
      kga = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "h") == 0) {
      h = force->numeric(FLERR, arg[iarg+1]);
      if (h <= 0)
	error->all(FLERR, "Henry's law solubility constant (H) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "temp") == 0) {
      temp = force->numeric(FLERR, arg[iarg+1]);
      if (temp <= 0)
	error->all(FLERR, "Temperature (temp) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor_vhead") == 0) {
      reactor_vhead = force->numeric(FLERR, arg[iarg+1]);
      if (reactor_vhead <= 0)
	error->all(FLERR, "Reactor headspace volume (reactor_vhead) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "reactor_pres") == 0) {
      reactor_pres = force->numeric(FLERR, arg[iarg+1]);
      if (reactor_pres <= 0)
	error->all(FLERR, "Reactor headspace pressure (reactor_pres) must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rg") == 0) {
      rg = force->numeric(FLERR, arg[iarg+1]);
      if (rg <= 0)
	error->all(FLERR, "Ideal gas constant must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "mw") == 0) {
      mw = force->numeric(FLERR, arg[iarg+1]);
      if (mw <= 0)
	error->all(FLERR, "Molar mass must be positive");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/gas_liquid command");
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixReactorGasLiquid::compute_scalar()
{
  double result = 0.0;

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      double p = grid->reac[igas][i];
      result += p;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);
  result *= rg * temp / reactor_pres;
  return result;
}


/* ---------------------------------------------------------------------- */

void FixReactorGasLiquid::compute()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double p_g2l, n_l2g;

  double vol = grid->cell_size * grid->cell_size * grid->cell_size;

  for (int i = 0; i < grid->ncells; i++) {
    if (!(grid->mask[i] & GHOST_MASK)) {
      p_g2l = kga * (grid->bulk[igas] - conc[iliquid][i]/(h * mw));
      n_l2g = -p_g2l * reactor_vhead / rg * temp;
      // update reaction rates
      reac[igas][i] += p_g2l;
      reac[iliquid][i] += n_l2g * mw / vol;
    }
  }
}