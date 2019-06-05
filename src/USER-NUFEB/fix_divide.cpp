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
#include "fix_divide.h"
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "force.h"
#include "lmptype.h"
#include "math_const.h"
#include "random_park.h"
#include "update.h"
#include "modify.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define DELTA 1.005

/* ---------------------------------------------------------------------- */

FixDivide::FixDivide(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6)
    error->all(FLERR, "Illegal fix nufeb/divide command");

  diameter = force->numeric(FLERR, arg[3]);
  eps_density = force->numeric(FLERR, arg[4]);
  seed = force->inumeric(FLERR, arg[5]);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  force_reneighbor = 1;
}

/* ---------------------------------------------------------------------- */

FixDivide::~FixDivide()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixDivide::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDivide::post_integrate()
{
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      double density = atom->rmass[i] /
	(4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]);

      if (atom->radius[i] * 2 >= diameter) {
        double newX, newY, newZ;

        double splitF = 0.4 + (random->uniform() * 0.2);
        double parentMass = atom->rmass[i] * splitF;
        double childMass = atom->rmass[i] - parentMass;

        double parentOuterMass = atom->outer_mass[i] * splitF;
        double childOuterMass = atom->outer_mass[i] - parentOuterMass;

        double thetaD = random->uniform() * 2 * MY_PI;
        double phiD = random->uniform() * (MY_PI);

        double oldX = atom->x[i][0];
        double oldY = atom->x[i][1];
        double oldZ = atom->x[i][2];

        // update parent
        atom->rmass[i] = parentMass;
        atom->outer_mass[i] = parentOuterMass;
        atom->radius[i] = pow(((6 * atom->rmass[i]) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
        atom->outer_radius[i] = pow((3.0 / (4.0 * MY_PI)) * ((atom->rmass[i] / density) + (parentOuterMass / eps_density)), (1.0 / 3.0));
        newX = oldX + (atom->outer_radius[i] * cos(thetaD) * sin(phiD) * DELTA);
        newY = oldY + (atom->outer_radius[i] * sin(thetaD) * sin(phiD) * DELTA);
        newZ = oldZ + (atom->outer_radius[i] * cos(phiD) * DELTA);
        if (newX - atom->outer_radius[i] < domain->boxlo[0]) {
          newX = domain->boxlo[0] + atom->outer_radius[i];
        } else if (newX + atom->outer_radius[i] > domain->boxhi[0]) {
          newX = domain->boxhi[0] - atom->outer_radius[i];
        }
        if (newY - atom->outer_radius[i] < domain->boxlo[1]) {
          newY = domain->boxlo[1] + atom->outer_radius[i];
        } else if (newY + atom->outer_radius[i] > domain->boxhi[1]) {
          newY = domain->boxhi[1] - atom->outer_radius[i];
        }
        if (newZ - atom->outer_radius[i] < domain->boxlo[2]) {
          newZ = domain->boxlo[2] + atom->outer_radius[i];
        } else if (newZ + atom->outer_radius[i] > domain->boxhi[2]) {
          newZ = domain->boxhi[2] - atom->outer_radius[i];
        }
        atom->x[i][0] = newX;
        atom->x[i][1] = newY;
        atom->x[i][2] = newZ;

        // create child
        double childRadius = pow(((6 * childMass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
        double childOuterRadius = pow((3.0 / (4.0 * MY_PI)) * ((childMass / density) + (childOuterMass / eps_density)), (1.0 / 3.0));
        double* coord = new double[3];
        newX = oldX - (childOuterRadius * cos(thetaD) * sin(phiD) * DELTA);
        newY = oldY - (childOuterRadius * sin(thetaD) * sin(phiD) * DELTA);
        newZ = oldZ - (childOuterRadius * cos(phiD) * DELTA);
        if (newX - childOuterRadius < domain->boxlo[0]) {
          newX = domain->boxlo[0] + childOuterRadius;
        } else if (newX + childOuterRadius > domain->boxhi[0]) {
          newX = domain->boxhi[0] - childOuterRadius;
        }
        if (newY - childOuterRadius < domain->boxlo[1]) {
          newY = domain->boxlo[1] + childOuterRadius;
        } else if (newY + childOuterRadius > domain->boxhi[1]) {
          newY = domain->boxhi[1] - childOuterRadius;
        }
        if (newZ - childOuterRadius < domain->boxlo[2]) {
          newZ = domain->boxlo[2] + childOuterRadius;
        } else if (newZ + childOuterRadius > domain->boxhi[2]) {
          newZ = domain->boxhi[2] - childOuterRadius;
        }
        coord[0] = newX;
        coord[1] = newY;
        coord[2] = newZ;

        int n = 0;
        atom->avec->create_atom(atom->type[i], coord);
        n = atom->nlocal - 1;

        atom->tag[n] = 0;
        atom->mask[n] = atom->mask[i];
        atom->image[n] = atom->image[i];

        atom->v[n][0] = atom->v[i][0];
        atom->v[n][1] = atom->v[i][1];
        atom->v[n][2] = atom->v[i][2];

        atom->omega[n][0] = atom->omega[i][0];
        atom->omega[n][1] = atom->omega[i][1];
        atom->omega[n][2] = atom->omega[i][2];

        atom->rmass[n] = childMass;
        atom->outer_mass[n] = childOuterMass;

        atom->radius[n] = childRadius;
        atom->outer_radius[n] = childOuterRadius;

        modify->create_attribute(n);

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

  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;
}
