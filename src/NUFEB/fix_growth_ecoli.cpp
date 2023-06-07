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

#include "fix_growth_ecoli.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include "atom.h"
#include "error.h"

#include "atom_vec_bacillus.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGrowthEcoli::FixGrowthEcoli(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 8)
    error->all(FLERR, "Illegal fix nufeb/growth/ecoli command");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/ecoli requires grid_style nufeb/chemostat");

  isuc = -1;
  io2 = -1;
  ico2 = -1;

  suc_affinity = 0.0;
  o2_affinity = 0.0;

  growth = 0.0;
  yield = 1.0;
  maintain = 0.0;
  decay = 0.0;

  avec = NULL;

  isuc = grid->find(arg[3]);
  if (isuc < 0)
    error->all(FLERR, "Can't find substrate name");
  suc_affinity = utils::numeric(FLERR,arg[4],true,lmp);
  if (suc_affinity <= 0)
    error->all(FLERR, "Sucrose affinity must be greater than zero");

  io2 = grid->find(arg[5]);
  if (io2 < 0)
    error->all(FLERR, "Can't find substrate(o2) name");
  o2_affinity = utils::numeric(FLERR,arg[6],true,lmp);
  if (o2_affinity <= 0)
    error->all(FLERR, "O2 affinity must be greater than zero");

  ico2 = grid->find(arg[7]);
  if (ico2 < 0)
    error->all(FLERR, "Can't find substrate(co2) name");

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "maintain") == 0) {
      maintain = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/ecoli command");
    }
  }

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
}

/* ---------------------------------------------------------------------- */

void FixGrowthEcoli::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      // cyanobacterial growth rate based on light(sub) and co2
      double tmp1 = growth * conc[isuc][i] / (suc_affinity + conc[isuc][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
      // sucrose export-induced growth reduction
      double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

      // nutrient utilization
      reac[isuc][i] -= 1 / yield * tmp1 * dens[igroup][i];
      reac[io2][i] -= 0.399 * (tmp1 + tmp2) * dens[igroup][i];
      reac[ico2][i] += 0.2 * (tmp1 + tmp2) * dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthEcoli::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
    // cyanobacterial growth rate based on light(sub) and co2
    double tmp1 = growth * conc[isuc][i] / (suc_affinity + conc[isuc][i]) * conc[io2][i] / (o2_affinity + conc[io2][i]);
    // sucrose export-induced growth reduction
    double tmp2 = maintain * conc[io2][i] / (o2_affinity + conc[io2][i]);

    grid->growth[igroup][i][0] = tmp1 - tmp2 - decay;
  }

  if (atom->coccus_flag) {
    update_atoms_coccus();
  } else {
    update_atoms_bacillus(avec);
  }
}
