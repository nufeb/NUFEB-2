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
#include "fix_monod_aob.h"
#include "atom.h"
#include "force.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixMonodAOB::FixMonodAOB(LAMMPS *lmp, int narg, char **arg) :
  FixMonod(lmp, narg, arg)
{
  if (narg < 8)
    error->all(FLERR, "Illegal fix nufeb/monod/aob command");

  dynamic_group_allow = 1;

  inh4 = -1;
  io2 = -1;
  ino2 = -1;

  nh4_affinity = 0.0;
  o2_affinity = 0.0;

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;

  inh4 = grid->find(arg[3]);
  if (inh4 < 0)
    error->all(FLERR, "Can't find substrate name");
  nh4_affinity = force->numeric(FLERR, arg[4]);

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate name");
  o2_affinity = force->numeric(FLERR, arg[6]);

  ino2 = grid->find(arg[7]);
  if (ino2 < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maintain") == 0) {
      maintain = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/monod/aob command");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonodAOB::compute()
{
  if (reaction_flag && growth_flag) {
    update_cells<1, 1>();
    update_atoms();
  } else if (reaction_flag && !growth_flag) {
    update_cells<1, 0>();
  } else if (!reaction_flag && growth_flag) {
    update_cells<0, 1>();
    update_atoms();
  }
}

/* ---------------------------------------------------------------------- */

template <int Reaction, int Growth>
void FixMonodAOB::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    double tmp1 = growth * conc[inh4][i] / (nh4_affinity + conc[inh4][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
    double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

    if (Reaction &&  !(grid->mask[i] & GHOST_MASK)) {
      reac[inh4][i] -= 1 / yield * tmp1 * dens[igroup][i];
      reac[io2][i] -= (4.57 - yield) / yield * tmp1 * dens[igroup][i] + tmp2 * dens[igroup][i];
      reac[ino2][i] += 1 / yield * tmp1 * dens[igroup][i];
    }

    if (Growth) {
      grid->growth[igroup][i][0] = tmp1 - tmp2 - decay;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMonodAOB::update_atoms()
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *biomass = atom->biomass;
  double *outer_radius = atom->outer_radius;
  double *outer_mass = atom->outer_mass;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell = grid->cell(x[i]);
      const double density = rmass[i] /
	(four_thirds_pi * radius[i] * radius[i] * radius[i]);
      double growth = grid->growth[igroup][cell][0];
      double ratio = rmass[i] / biomass[i];
      // forward Eular to update biomass and rmass
      biomass[i] = biomass[i] * (1 + growth * dt);
      rmass[i] = rmass[i] * (1 + growth * dt * ratio);
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      outer_mass[i] = 0;
      outer_radius[i] = radius[i];
    }
  }
}
