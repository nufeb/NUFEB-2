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
#include "atom.h"
#include "error.h"

#include "atom_vec_bacillus.h"
#include "fix_growth_simple.h"
#include "grid.h"
#include "group.h"
#include "modify.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixGrowthSimple::FixGrowthSimple(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/growth/monod command");

  dynamic_group_allow = 1;

  isub = -1;

  growth = 0.0;

  avec = nullptr;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate name");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/simple command");
    }
  }
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
}

/* ---------------------------------------------------------------------- */

void FixGrowthSimple::compute()
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
void FixGrowthSimple::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (Growth) {
      grid->growth[igroup][i][0] = growth;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthSimple::update_atoms()
{
  if (atom->coccus_flag) {
    update_atoms_coccus();
  } else {
    update_atoms_bacillus(avec);
  }
}
