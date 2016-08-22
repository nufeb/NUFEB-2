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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom_vec_bio.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 500000

/* ---------------------------------------------------------------------- */

AtomVecBio::AtomVecBio(LAMMPS *lmp) : AtomVecSphere(lmp)
{
  size_data_atom = 8;

  atom->outerMass = NULL;
  atom->outerRadius = NULL;
  atom->virtualMass = NULL;
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::init()
{
  AtomVecSphere::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"growth") == 0 || strcmp(modify->fix[i]->style,"divide") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      radvary = 1;
      comm_x_only = 0;
      size_forward = 5;
    }

  int ntypes = atom->ntypes;
  virtualMass = memory->grow(atom->virtualMass,ntypes+1,"atom:virtualMass");
  for (int i = 0; i < ntypes+1; i++) {
    virtualMass[i] = 0;
  }
}

void AtomVecBio::grow(int n)
{
  AtomVecSphere::grow(n);
  outerMass = memory->grow(atom->outerMass,nmax,"atom:outerMass");
  outerRadius = memory->grow(atom->outerRadius,nmax,"atom:outerRadius");
}

void AtomVecBio::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  AtomVecSphere::create_atom(itype, coord);

  outerRadius[nlocal] = 0.5;
  outerMass[nlocal] = 0.0;
}


void AtomVecBio::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  AtomVecSphere::data_atom(coord, imagetmp, values);

  outerRadius[nlocal] = 0.5 * atof(values[7]);
  if (type[nlocal] != 1 && outerRadius[nlocal] != radius[nlocal]) {
    error->one(FLERR,"Outer radius must be equal to radius for all AOB, NOB, EPS, and inert particles");
  }
  else if (outerRadius[nlocal] < radius[nlocal]) {
    error->one(FLERR,"Outer radius must be greater than or equal to radius");
  }

  outerMass[nlocal] = 0.0;

  //outerMass[nlocal] = (4.0*MY_PI/3.0)*((outerRadius[nlocal]*outerRadius[nlocal]*outerRadius[nlocal])-(radius[nlocal]*radius[nlocal]*radius[nlocal]))*30;
}

void AtomVecBio::grow_reset()
{
  AtomVecSphere::grow_reset();

  outerMass = atom->outerMass;
  outerRadius = atom->outerRadius;
}

void AtomVecBio::copy(int i, int j, int delflag)
{
  outerMass[j] = outerMass[i];
  outerRadius[j] = outerRadius[i];

  AtomVecSphere::copy(i, j, delflag);

}

bigint AtomVecBio::memory_usage()
{
  bigint bytes = AtomVecSphere::memory_usage();

  if (atom->memcheck("outerMass")) bytes += memory->usage(outerMass,nmax);
  if (atom->memcheck("outerRadius")) bytes += memory->usage(outerRadius,nmax);

  return bytes;
}
