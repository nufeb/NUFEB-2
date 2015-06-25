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

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecBio::AtomVecBio(LAMMPS *lmp) : AtomVecSphere(lmp)
{
  size_data_atom = 13;
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

  virtualMass = memory->grow(atom->virtualMass,3,"atom:virtualMass");
  for (int i = 0; i < 3; i++) {
    virtualMass[i] = 0;
  }
}

void AtomVecBio::grow(int n)
{
  AtomVecSphere::grow(n);
  sub = memory->grow(atom->sub,nmax,"atom:sub");
  o2 = memory->grow(atom->o2,nmax,"atom:o2");
  nh4 = memory->grow(atom->nh4,nmax,"atom:nh4");
  no2 = memory->grow(atom->no2,nmax,"atom:no2");
  no3 = memory->grow(atom->no3,nmax,"atom:no3");
  outerRadius = memory->grow(atom->outerRadius,nmax,"atom:outerRadius");
  outerMass = memory->grow(atom->outerMass,nmax,"atom:outerMass");
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
  sub[nlocal] = atof(values[8]);
  o2[nlocal] = atof(values[9]);
  nh4[nlocal] = atof(values[10]);
  no2[nlocal] = atof(values[11]);
  no3[nlocal] = atof(values[12]);

  outerMass[nlocal] = 0.0;

  //outerMass[nlocal] = (4.0*MY_PI/3.0)*((outerRadius[nlocal]*outerRadius[nlocal]*outerRadius[nlocal])-(radius[nlocal]*radius[nlocal]*radius[nlocal]))*30;
}

void AtomVecBio::copy(int i, int j, int delflag)
{
  sub[j] = sub[i];
  o2[j] = o2[i];
  nh4[j] = nh4[i];
  no2[j] = no2[i];
  no3[j] = no3[i];
  outerMass[j] = outerMass[i];
  outerRadius[j] = outerRadius[i];

  AtomVecSphere::copy(i, j, delflag);

}
