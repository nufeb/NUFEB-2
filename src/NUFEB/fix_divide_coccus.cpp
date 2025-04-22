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

#include "fix_divide_coccus.h"

#include "atom.h"
#include "atom_masks.h"
#include "atom_vec.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_property_cycletime.h"
#include "grid.h"
#include "lmptype.h"
#include "math_const.h"
#include "modify.h"
#include "random_park.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define DELTA 1.005

/* ---------------------------------------------------------------------- */

FixDivideCoccus::FixDivideCoccus(LAMMPS *lmp, int narg, char **arg) : FixDivide(lmp, narg, arg)
{
  if (!atom->coccus_flag) error->all(FLERR, "fix nufeb/division/coccus requires coccus atom style");

  if (narg < 6) error->all(FLERR, "Illegal fix nufeb/division/coccus command");

  max_dia = 1.0;
  max_time = 1.0;
  eps_density = 30;
  x_flag = y_flag = z_flag = 1;
  size_flag = time_flag = 0;
  fix_flag = 0;
  fix_ct = nullptr;

  group_id = new char[strlen(arg[1])+1];
  strcpy(group_id, arg[1]);

  if (strcmp(arg[3], "size") == 0) {
    size_flag = 1;
    max_dia = utils::numeric(FLERR, arg[4], true, lmp);
  } else if (strcmp(arg[3], "time") == 0) {
    time_flag = 1;
    max_time = utils::numeric(FLERR, arg[4], true, lmp);
  } else{
    error->all(FLERR, "Illegal fix nufeb/division/coccus command");
  }

  seed = utils::inumeric(FLERR, arg[5], true, lmp);

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "epsdens") == 0) {
      eps_density = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    }  else if (strcmp(arg[iarg], "x") == 0) {
      x_flag = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "y") == 0) {
      y_flag = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "z") == 0) {
      z_flag = utils::inumeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "fix") == 0) {
      fix_flag = 1;
      iarg += 1;
    } else {
      error->all(FLERR, "Illegal fix nufeb/division/coccus command");
    }
  }
  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
}

/* ---------------------------------------------------------------------- */

FixDivideCoccus::~FixDivideCoccus()
{
  delete random;
}

/* ---------------------------------------------------------------------- */
void FixDivideCoccus::post_constructor()
{
  if (time_flag) {
    // create fix nufeb/property/cycletime
    int nargs = 5;
    char **fixarg = new char *[nargs];
    std::string max_time_str = std::to_string(max_time);
    fixarg[0] = (char *) "div_coccus";
    fixarg[1] = group_id;
    fixarg[2] = (char *) "nufeb/property/cycletime";
    fixarg[3] = (char *) "max_time";
    fixarg[4] = strdup(max_time_str.c_str());
    ;

    modify->add_fix(nargs, fixarg, 1);
    delete[] fixarg;
    fix_ct = (FixPropertyCycletime *) modify->fix[modify->nfix - 1];
  }
}
/* ---------------------------------------------------------------------- */

void FixDivideCoccus::compute()
{
  double **x = atom->x;
  int nlocal = atom->nlocal;
  const double third = 1.0 / 3.0;
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int divide = 0;
      if (size_flag) {
        if (atom->radius[i] * 2 > max_dia) divide = 1;
      } else if (time_flag) {
        // convert growth rate to doubling time
        const int cell = grid->cell(x[i]);
        double growth = grid->growth[igroup][cell][0];
        double div_time = log(2) / growth;
        if (fix_ct->aprop[i][0] > div_time) divide = 1;
      }

      if (divide) {
        double imass, jmass, iouter_mass, jouter_mass, density;

        density = atom->rmass[i] /
            (four_thirds_pi * atom->radius[i] * atom->radius[i] * atom->radius[i]);

        if (fix_flag) {
          imass = atom->rmass[i];
          jmass = atom->rmass[i];

          iouter_mass = atom->outer_mass[i];
          jouter_mass = atom->outer_mass[i];
        } else {
          double split = 0.4 + (random->uniform() * 0.2);
          imass = atom->rmass[i] * split;
          jmass = atom->rmass[i] - imass;

          iouter_mass = atom->outer_mass[i] * split;
          jouter_mass = atom->outer_mass[i] - iouter_mass;
        }

        double theta = random->uniform() * 2 * MY_PI;
        double phi = random->uniform() * (MY_PI);

        double oldx = atom->x[i][0];
        double oldy = atom->x[i][1];
        double oldz = atom->x[i][2];

        // update daughter cell i
        atom->rmass[i] = imass;
        atom->outer_mass[i] = iouter_mass;
        atom->radius[i] = pow(((6 * atom->rmass[i]) / (density * MY_PI)), third) * 0.5;
        atom->outer_radius[i] = pow(
            three_quarters_pi * ((atom->rmass[i] / density) + (iouter_mass / eps_density)), third);

        if (x_flag) {
          double newx = oldx + (atom->radius[i] * cos(theta) * sin(phi) * DELTA);
          if (newx - atom->radius[i] < domain->boxlo[0])
            newx = domain->boxlo[0] + atom->radius[i];
          else if (newx + atom->radius[i] > domain->boxhi[0])
            newx = domain->boxhi[0] - atom->radius[i];
          atom->x[i][0] = newx;
        }
        if (y_flag) {
          double newy = oldy + (atom->radius[i] * sin(theta) * sin(phi) * DELTA);
          if (newy - atom->radius[i] < domain->boxlo[1])
            newy = domain->boxlo[1] + atom->radius[i];
          else if (newy + atom->radius[i] > domain->boxhi[1])
            newy = domain->boxhi[1] - atom->radius[i];
          atom->x[i][1] = newy;
        }
        if (z_flag) {
          double newz = oldz + (atom->radius[i] * cos(phi) * DELTA);
          if (newz - atom->radius[i] < domain->boxlo[2])
            newz = domain->boxlo[2] + atom->radius[i];
          else if (newz + atom->radius[i] > domain->boxhi[2])
            newz = domain->boxhi[2] - atom->radius[i];
          atom->x[i][2] = newz;
        }

        // create daughter cell j
        double jradius = pow(((6 * jmass) / (density * MY_PI)), third) * 0.5;
        double jouter_radius =
            pow(three_quarters_pi * ((jmass / density) + (jouter_mass / eps_density)), third);
        double *coord = new double[3];

        if (x_flag) {
          double newx = oldx - (jradius * cos(theta) * sin(phi) * DELTA);
          if (newx - jradius < domain->boxlo[0])
            newx = domain->boxlo[0] + jradius;
          else if (newx + jradius > domain->boxhi[0])
            newx = domain->boxhi[0] - jradius;
          coord[0] = newx;
        } else {
          coord[0] = atom->x[i][0];
        }

        if (y_flag) {
          double newy = oldy - (jradius * sin(theta) * sin(phi) * DELTA);
          if (newy - jradius < domain->boxlo[1])
            newy = domain->boxlo[1] + jouter_radius;
          else if (newy + jradius > domain->boxhi[1])
            newy = domain->boxhi[1] - jradius;
          coord[1] = newy;
        } else {
          coord[1] = atom->x[i][1];
        }

        if (z_flag) {
          double newz = oldz - (jradius * cos(phi) * DELTA);
          if (newz - jradius < domain->boxlo[2])
            newz = domain->boxlo[2] + jradius;
          else if (newz + jradius > domain->boxhi[2])
            newz = domain->boxhi[2] - jradius;
          coord[2] = newz;
        } else {
          coord[2] = atom->x[i][2];
        }

        atom->avec->create_atom(atom->type[i], coord);
        int j = atom->nlocal - 1;

        atom->tag[j] = 0;
        atom->mask[j] = atom->mask[i];
        atom->v[j][0] = atom->v[i][0];
        atom->v[j][1] = atom->v[i][1];
        atom->v[j][2] = atom->v[i][2];
        atom->f[j][0] = atom->f[i][0];
        atom->f[j][1] = atom->f[i][1];
        atom->f[j][2] = atom->f[i][2];
        atom->omega[j][0] = atom->omega[i][0];
        atom->omega[j][1] = atom->omega[i][1];
        atom->omega[j][2] = atom->omega[i][2];
        atom->rmass[j] = jmass;
        atom->biomass[j] = atom->biomass[i];
        atom->radius[j] = jradius;
        atom->outer_mass[j] = jouter_mass;
        atom->outer_radius[j] = jouter_radius;

        modify->create_attribute(j);

        for (int m = 0; m < modify->nfix; m++) modify->fix[m]->update_arrays(i, j);

        delete[] coord;
      }
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT) error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable) atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}
