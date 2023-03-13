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
#include "update.h"
#include "atom.h"
#include "grid.h"
#include "math_const.h"
#include "atom_vec_bacillus.h"
#include "fix_growth.h"
#include "grid_masks.h"
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define MIN_MASS 1e-30

/* ---------------------------------------------------------------------- */

FixGrowth::FixGrowth(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  dt = 1.0;
  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

int FixGrowth::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "nevery") == 0) {
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery <= 0) error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixGrowth::init()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixGrowth::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

int FixGrowth::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  mask |= CHEMISTRY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGrowth::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  update_atoms();
}

/* ---------------------------------------------------------------------- */

void FixGrowth::chemistry_nufeb()
{
  update_cells();
}

/* ---------------------------------------------------------------------- */

void FixGrowth::update_atoms_coccus()
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_radius = atom->outer_radius;
  double *outer_mass = atom->outer_mass;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;
  int mass_flag = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      const int cell = grid->cell(x[i]);
      // skip atoms in ghost cells
      if (grid->mask[cell] & GHOST_MASK) continue;

      const double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);
      double growth = grid->growth[igroup][cell][0];
      // forward Euler to update rmass
      rmass[i] = rmass[i] * (1 + growth * dt);
      if (rmass[i] <= 0) {
        rmass[i] = MIN_MASS;
        mass_flag = 1;
      }

      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
      outer_mass[i] = 0;
      outer_radius[i] = radius[i];
    }
  }

  if (mass_flag)
    error->warning(FLERR,"Negative atom mass, reset value to 1e-30kg. Consider using fix nufeb/death/diameter");
}

/* ---------------------------------------------------------------------- */

void FixGrowth::update_atoms_bacillus(AtomVecBacillus *&avec)
{
  double **x = atom->x;
  double *rmass = atom->rmass;

  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  int mass_flag = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
      double acircle = MY_PI*atom->radius[i]*atom->radius[i];

      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];
      double length = bonus->length;

      double new_length;
      const int cell = grid->cell(x[i]);
      // skip atoms in ghost cells
      if (grid->mask[cell] & GHOST_MASK) continue;

      const double density = rmass[i] /	(vsphere + acircle * bonus->length);
      double growth = grid->growth[igroup][cell][0];
      // forward Eular to update rmass
      rmass[i] = rmass[i] * (1 + growth * dt);
      if (rmass[i] <= 0) {
        rmass[i] = MIN_MASS;
        mass_flag = 1;
      }

      new_length = (rmass[i] - density * vsphere) / (density * acircle);
      bonus->length = new_length;
      // update coordinates of two poles
      double *pole1 = bonus->pole1;
      double *pole2 = bonus->pole2;

      pole1[0] *= new_length/length;
      pole1[1] *= new_length/length;
      pole1[2] *= new_length/length;
      pole2[0] *= new_length/length;
      pole2[1] *= new_length/length;
      pole2[2] *= new_length/length;
    }
  }

  if (mass_flag)
    error->warning(FLERR,"Negative atom mass, reset value to 1e-30kg. Consider using fix nufeb/death/diameter");
}
