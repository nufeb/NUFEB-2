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

#include "fix_bio_eps_extract.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "atom_vec_bio.h"
#include "domain.h"
#include "error.h"

#include "fix_bio_fluid.h"
#include "force.h"
#include "input.h"
#include "lmptype.h"
#include "math_const.h"
#include "pointers.h"
#include "random_park.h"
#include "update.h"
#include "variable.h"
#include "group.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001
#define DELTA 1.005

// enum{PAIR,KSPACE,ATOM};
// enum{DIAMETER,CHARGE};

/* ---------------------------------------------------------------------- */

FixEPSExtract::FixEPSExtract(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (narg < 7)
    error->all(FLERR, "Illegal fix eps_extract command: not enough arguments");

  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix eps_extract command: nevery is negative");

  var = new char*[2];
  ivar = new int[2];

  int i;
  for (i = 0; i < 2; i++) {
    int n = strlen(&arg[4 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[4 + i][2]);
  }

  seed = atoi(arg[6]);
  if (seed <= 0)
    error->all(FLERR, "Illegal fix eps extract command: seed should be greater than 0");

  demflag = 0;
  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (demflag != 0 && demflag != 1)
        error->all(FLERR, "Illegal fix eps_extract command: demflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix eps_extract command");
  }

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  // Set up renieghbouring here: required for re building the neighbour list: fix pour/ deposit

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

FixEPSExtract::~FixEPSExtract() {
  delete random;
  int i;
  for (i = 0; i < 2; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;
}

/* ---------------------------------------------------------------------- */

int FixEPSExtract::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixEPSExtract::init() {
  // fprintf(stdout, "called once?\n");
  if (!atom->radius_flag)
    error->all(FLERR, "Fix eps extract requires atom attribute diameter");

  int i;
  for (i = 0; i < 2; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR, "Variable name for fix eps extract does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR, "Variable for fix eps extract is invalid style");
  }

  eps_ratio = input->variable->compute_equal(ivar[0]);
  eps_density = input->variable->compute_equal(ivar[1]);

  int id = 0;
  for (int i = 1; i < group->ngroup; i++) {
    if (strcmp(group->names[i], "HET") == 0) {
      id = i;
      break;
    }
  }

  if (groupbit != pow(2, id)) {
    error->all(FLERR, "Fix eps extract is only valid for particles of type HET");
  }

  eps_type = avec->eps_type;

  if (eps_type == 0) {
    error->all(FLERR, "Cannot find EPS type");
  }

  nufebFoam = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "nufebFoam") == 0) {
      nufebFoam = static_cast<FixFluid *>(lmp->modify->fix[j]);
      break;
    }
  }
}

void FixEPSExtract::post_integrate() {
  if (nevery == 0)
    return;
  if (update->ntimestep % nevery)
    return;
  if (nufebFoam != NULL && nufebFoam->demflag)
    return;
  if (demflag)
    return;

  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {

      if ((avec->outer_radius[i] / atom->radius[i]) > eps_ratio) {
        avec->outer_mass[i] = (4.0 * MY_PI / 3.0) * ((avec->outer_radius[i] * avec->outer_radius[i] * avec->outer_radius[i]) - (atom->radius[i] * atom->radius[i] * atom->radius[i])) * eps_density;

        double splitF = 0.4 + (random->uniform() * 0.2);

        double newOuterMass = avec->outer_mass[i] * splitF;
        double eps_mass = avec->outer_mass[i] - newOuterMass;

        avec->outer_mass[i] = newOuterMass;

        double density = atom->rmass[i] / (4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]);
        avec->outer_radius[i] = pow((3.0 / (4.0 * MY_PI)) * ((atom->rmass[i] / density) + (avec->outer_mass[i] / eps_density)), (1.0 / 3.0));

        double thetaD = random->uniform() * 2 * MY_PI;
        double phiD = random->uniform() * (MY_PI);

        double oldX = atom->x[i][0];
        double oldY = atom->x[i][1];
        double oldZ = atom->x[i][2];

        //create child
        double childRadius = pow(((6 * eps_mass) / (eps_density * MY_PI)), (1.0 / 3.0)) * 0.5;
        double* coord = new double[3];
        double newX = oldX - ((childRadius + avec->outer_radius[i]) * cos(thetaD) * sin(phiD) * DELTA);
        double newY = oldY - ((childRadius + avec->outer_radius[i]) * sin(thetaD) * sin(phiD) * DELTA);
        double newZ = oldZ - ((childRadius + avec->outer_radius[i]) * cos(phiD) * DELTA);
        if (newX - childRadius < xlo) {
          newX = xlo + childRadius;
        } else if (newX + childRadius > xhi) {
          newX = xhi - childRadius;
        }
        if (newY - childRadius < ylo) {
          newY = ylo + childRadius;
        } else if (newY + childRadius > yhi) {
          newY = yhi - childRadius;
        }
        if (newZ - childRadius < zlo) {
          newZ = zlo + childRadius;
        } else if (newZ + childRadius > zhi) {
          newZ = zhi - childRadius;
        }
        coord[0] = newX;
        coord[1] = newY;
        coord[2] = newZ;

        int n = 0;
        atom->avec->create_atom(eps_type, coord);
        n = atom->nlocal - 1;

        atom->tag[n] = 0;
        atom->mask[n] = avec->eps_mask;

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
        avec->outer_mass[n] = 0;

        atom->torque[n][0] = atom->torque[i][0];
        atom->torque[n][1] = atom->torque[i][1];
        atom->torque[n][2] = atom->torque[i][2];
        atom->radius[n] = childRadius;
        avec->outer_radius[n] = childRadius;

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

/* ---------------------------------------------------------------------- */

int FixEPSExtract::modify_param(int narg, char **arg) {
  if (strcmp(arg[0], "demflag") == 0) {
    if (narg != 2)
      error->all(FLERR, "Illegal fix_modify command");
    demflag = force->inumeric(FLERR, arg[1]);
    if (demflag != 1 && demflag != 0)
      error->all(FLERR, "Illegal fix_modify command: demflag");
    return 2;
  }
  return 0;
}

