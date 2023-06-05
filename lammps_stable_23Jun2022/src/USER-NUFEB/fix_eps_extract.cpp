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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_eps_extract.h"
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "math_const.h"
#include "random_park.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "group.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define DELTA 1.005

/* ---------------------------------------------------------------------- */

FixEPSExtract::FixEPSExtract(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8)
    error->all(FLERR, "Illegal fix nufeb/eps_extract command");

  type = utils::inumeric(FLERR,arg[3],true,lmp);
  ieps = group->find(arg[4]);
  if (ieps < 0)
    error->all(FLERR, "Can't find group in fix nufeb/eps_extract command");
  ratio = utils::numeric(FLERR,arg[5],true,lmp);
  eps_density = utils::numeric(FLERR,arg[6],true,lmp);
  seed = utils::inumeric(FLERR,arg[7],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  force_reneighbor = 1;
}

/* ---------------------------------------------------------------------- */

FixEPSExtract::~FixEPSExtract()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixEPSExtract::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixEPSExtract::modify_param(int narg, char **arg)
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

void FixEPSExtract::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixEPSExtract::post_neighbor()
{
  // reset reneighbor flag
  next_reneighbor = 0;
}

/* ---------------------------------------------------------------------- */

void FixEPSExtract::compute()
{
  int nlocal = atom->nlocal;
  int eps_mask = group->bitmask[ieps];
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if ((atom->outer_radius[i] / atom->radius[i]) > ratio) {
        atom->outer_mass[i] = (4.0 * MY_PI / 3.0) *
            ((atom->outer_radius[i] * atom->outer_radius[i] * atom->outer_radius[i]) -
        	(atom->radius[i] * atom->radius[i] * atom->radius[i])) * eps_density;

        double split = 0.4 + (random->uniform() * 0.2);

        double new_outer_mass = atom->outer_mass[i] * split;
        double eps_mass = atom->outer_mass[i] - new_outer_mass;

        atom->outer_mass[i] = new_outer_mass;

        double density = atom->rmass[i] / (4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]);
        atom->outer_radius[i] = pow((3.0 / (4.0 * MY_PI)) * ((atom->rmass[i] / density) + (atom->outer_mass[i] / eps_density)), (1.0 / 3.0));

        double theta = random->uniform() * 2 * MY_PI;
        double phi = random->uniform() * (MY_PI);

        double oldx = atom->x[i][0];
        double oldy = atom->x[i][1];
        double oldz = atom->x[i][2];

        // create EPS atom
        double eps_radius = pow(((6 * eps_mass) / (eps_density * MY_PI)), (1.0 / 3.0)) * 0.5;
        double *coord = new double[3];
        double newx = oldx - ((eps_radius + atom->outer_radius[i]) * cos(theta) * sin(phi) * DELTA);
        double newy = oldy - ((eps_radius + atom->outer_radius[i]) * sin(theta) * sin(phi) * DELTA);
        double newz = oldz - ((eps_radius + atom->outer_radius[i]) * cos(phi) * DELTA);
        if (newx - eps_radius < domain->boxlo[0]) {
          newx = domain->boxlo[0] + eps_radius;
        } else if (newx + eps_radius > domain->boxhi[0]) {
          newx = domain->boxhi[0] - eps_radius;
        }
        if (newy - eps_radius < domain->boxlo[1]) {
          newy = domain->boxlo[1] + eps_radius;
        } else if (newy + eps_radius > domain->boxhi[1]) {
          newy = domain->boxhi[1] - eps_radius;
        }
        if (newz - eps_radius < domain->boxlo[2]) {
          newz = domain->boxlo[2] + eps_radius;
        } else if (newz + eps_radius > domain->boxhi[2]) {
          newz = domain->boxhi[2] - eps_radius;
        }
        coord[0] = newx;
        coord[1] = newy;
        coord[2] = newz;

        atom->avec->create_atom(type, coord);
        int n = atom->nlocal - 1;

        atom->tag[n] = 0;
        atom->mask[n] = 1 | eps_mask;
        atom->v[n][0] = atom->v[i][0];
        atom->v[n][1] = atom->v[i][1];
        atom->v[n][2] = atom->v[i][2];
        atom->f[n][0] = atom->f[i][0];
        atom->f[n][1] = atom->f[i][1];
        atom->f[n][2] = atom->f[i][2];
        atom->omega[n][0] = atom->omega[i][0];
        atom->omega[n][1] = atom->omega[i][1];
        atom->omega[n][2] = atom->omega[i][2];
        atom->rmass[n] = eps_mass;
        atom->biomass[n] = eps_mass;
        atom->radius[n] = eps_radius;
        atom->outer_mass[n] = 0.0;
        atom->outer_radius[n] = eps_radius;

        modify->create_attribute(n);

        for (int m = 0; m < modify->nfix; m++)
          modify->fix[m]->update_arrays(i, n);

        delete[] coord;
      }
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}
