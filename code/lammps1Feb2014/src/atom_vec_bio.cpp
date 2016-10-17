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
  size_data_atom = 9;

  //atom
  atom->outerMass = NULL;
  atom->outerRadius = NULL;
  atom->virtualMass = NULL;
  atom->atom_growth = NULL;

  //type
  atom->ks = NULL;
  atom->growth = NULL;
  atom->yield = NULL;
  atom->typeName = NULL;
  atom->anabCoeff = NULL;
  atom->catCoeff = NULL;

  //nutrient
  atom->nNutrients = 0;
  atom->iniS = NULL;
  atom->nuName = NULL;
  atom->diffCoeff = NULL;
  atom->nuS = NULL;
  atom->nuR = NULL;
  atom->nuType = NULL;

}

/* ---------------------------------------------------------------------- */

void AtomVecBio::init()
{
  AtomVecSphere::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"growth") == 0 || strcmp(modify->fix[i]->style,"divide") == 0 || strcmp(modify->fix[i]->style,"metabolism") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      radvary = 1;
      comm_x_only = 0;
      size_forward = 5;
    }
}

void AtomVecBio::grow(int n)
{
  AtomVecSphere::grow(n);
  outerMass = memory->grow(atom->outerMass,nmax,"atom:outerMass");
  outerRadius = memory->grow(atom->outerRadius,nmax,"atom:outerRadius");
  atom_growth = memory->grow(atom->atom_growth,nmax,"atom:atom_growth");
}

void AtomVecBio::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  AtomVecSphere::create_atom(itype, coord);

  outerRadius[nlocal] = 0.0;
  outerMass[nlocal] = 0.0;
  atom_growth[nlocal] = 0.0;
}


void AtomVecBio::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  typeName = atom->typeName;

  AtomVecSphere::data_atom(coord, imagetmp, values);

  outerRadius[nlocal] = 0.5 * atof(values[7]);
  //outerRadius[nlocal] = 0.5 * atof(values[7]);
//  if (type[nlocal] != 1 && outerRadius[nlocal] != radius[nlocal]) {
//    error->one(FLERR,"Outer radius must be equal to radius for all AOB, NOB, EPS, and inert particles");
//  }
  if (outerRadius[nlocal] < radius[nlocal]) {
    error->one(FLERR,"Outer radius must be greater than or equal to radius");
  }

  outerMass[nlocal] = (4.0*MY_PI/3.0)*((outerRadius[nlocal]*outerRadius[nlocal]*outerRadius[nlocal])
      -(radius[nlocal]*radius[nlocal]*radius[nlocal]))*30;

  char *name;
  int type = atom->type[nlocal];

  int n = strlen(values[8]) + 1;
  name = new char[n];
  strcpy(name,values[8]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(name[i]) && name[i] != '_')
      error->all(FLERR,"Type name must be "
                 "alphanumeric or underscore characters");

  for (int i = 0; i < ntypes+1; i++)
    if (typeName[i] != NULL)
        if((strcmp(typeName[i], name) == 0) && (i != type)){
      error->one(FLERR,"Repeat type names");
    }

  if (typeName[type] == NULL){
    typeName[type] = new char[n];
  }
  else if (strcmp(typeName[type], name) != 0){
    error->one(FLERR,"Incompatible type names");
  }

  strcpy(typeName[type],name);

  delete[] name;
  //printf("name = %s type = %i \n", typeName[type], type);
}

void AtomVecBio::grow_reset()
{
  AtomVecSphere::grow_reset();

  outerMass = atom->outerMass;
  outerRadius = atom->outerRadius;
  atom_growth = atom->atom_growth;
}

void AtomVecBio::copy(int i, int j, int delflag)
{
  outerMass[j] = outerMass[i];
  outerRadius[j] = outerRadius[i];
  atom_growth[j] = atom_growth[i];

  AtomVecSphere::copy(i, j, delflag);

}

bigint AtomVecBio::memory_usage()
{
  bigint bytes = AtomVecSphere::memory_usage();

  if (atom->memcheck("outerMass")) bytes += memory->usage(outerMass,nmax);
  if (atom->memcheck("outerRadius")) bytes += memory->usage(outerRadius,nmax);
  if (atom->memcheck("atom_growth")) bytes += memory->usage(atom_growth,nmax);

  return bytes;
}

