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

#include "atom_vec_bio.h"

#include <stdlib.h>
#include <string.h>
#include <cctype>
#include <math.h>

#include "atom.h"
#include "error.h"
#include "fix_adapt.h"
#include "lmptype.h"
#include "math_const.h"
#include "memory.h"

#include "bio.h"
#include "modify.h"
#include "pointers.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define DELTA 500000

/* ---------------------------------------------------------------------- */

AtomVecBio::AtomVecBio(LAMMPS *lmp) : AtomVecSphere(lmp)
{
  size_data_atom = 8;

  //atom
  outerMass = memory->create(outerMass,nmax,"atom:outerMass");
  outerRadius = memory->create(outerRadius,nmax,"atom:outerRadius");;
  //atom_q = memory->create(atom_q,nmax,"atom:atom_q");
  typeEPS = 0;
  typeDEAD = 0;
  maskEPS = 0;
  maskHET = 0;
  maskDEAD = 0;
  //virtualMass = NULL;

  //instantiate BIO class
  bio = new BIO(lmp);
}

AtomVecBio::~AtomVecBio()
{
  delete bio;
  memory->destroy(outerMass);
  memory->destroy(outerRadius);
  //memory->destroy(atom_q);
  //memory->destroy(atom->virtualMass);
}

/* ---------------------------------------------------------------------- */

void AtomVecBio::init()
{
  AtomVecSphere::init();

  // set radvary if particle diameters are time-varying due to fix adapt

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"growth") == 0 || strcmp(modify->fix[i]->style,"divide") == 0 || strcmp(modify->fix[i]->style,"kinetics/monod") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      //need to fix here!
      //radvary = 1;
      comm_x_only = 0;
      size_forward = 5;
    }
  set_group_mask();
}

void AtomVecBio::grow(int n)
{
  AtomVecSphere::grow(n);
  outerMass = memory->grow(outerMass,nmax,"atom:outerMass");
  outerRadius = memory->grow(outerRadius,nmax,"atom:outerRadius");
  //atom_q = memory->grow(atom_q,nmax,"atom:atom_q");
}

void AtomVecBio::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  AtomVecSphere::create_atom(itype, coord);

  outerRadius[nlocal] = 0.0;
  outerMass[nlocal] = 0.0;
//  atom_q[nlocal] = 0.0;
}


void AtomVecBio::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  char **typeName  = bio->typeName;;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  AtomVecSphere::data_atom(coord, imagetmp, values);

  outerRadius[nlocal] = 0.5 * atof(values[7]);
  // outerRadius[nlocal] = 0.5 * atof(values[7]);
  // if (type[nlocal] != 1 && outerRadius[nlocal] != radius[nlocal]) {
  //   error->one(FLERR,"Outer radius must be equal to radius for all AOB, NOB, EPS, and inert particles");
  // }
  if (outerRadius[nlocal] < radius[nlocal]) {
    error->one(FLERR,"Outer radius must be greater than or equal to radius");
  }

  outerMass[nlocal] = (4.0*MY_PI/3.0)*((outerRadius[nlocal]*outerRadius[nlocal]*outerRadius[nlocal])
      -(radius[nlocal]*radius[nlocal]*radius[nlocal])) * 30;

  // double aa = 0.1367e-15;
  // printf("diameter = %e, \n", cbrt((aa*3)/(290*4*MY_PI)) * 2);

  // convert kg to mol
  // rmass[nlocal] = rmass[nlocal]/2.46e-2;
  // outerMass[nlocal] = outerMass[nlocal]/2.46e-2;

  // char *name;
  // int type = atom->type[nlocal];

  // int n = strlen(values[8]) + 1;
  // name = new char[n];
  // strcpy(name,values[8]);

  // for (int i = 0; i < n-1; i++)
  //   if (!isalnum(name[i]) && name[i] != '_')
  //     error->all(FLERR,"Type name must be "
  //                "alphanumeric or underscore characters");

  // for (int i = 0; i < ntypes+1; i++)
  //   if (typeName[i] != NULL)
  //       if((strcmp(typeName[i], name) == 0) && (i != type)){
  //     error->one(FLERR,"Repeat type names");
  //   }

  // if (typeName[type] == NULL){
  //   typeName[type] = new char[n];
  // }
  // else if (strcmp(typeName[type], name) != 0){
  //   error->one(FLERR,"Incompatible type names");
  // }

  // if (strcmp(name,"eps") == 0) {
  //   typeEPS = type;
  // } else if (strcmp(name,"dead") == 0) {
  //   typeDEAD = type;
  // }

  // strcpy(typeName[type],name);

  // delete[] name;
  // printf("name = %s type = %i \n", typeName[type], type);
}

void AtomVecBio::copy(int i, int j, int delflag)
{
  outerMass[j] = outerMass[i];
  outerRadius[j] = outerRadius[i];
  //atom_q[j] = atom_q[i];

  AtomVecSphere::copy(i, j, delflag);

}

bigint AtomVecBio::memory_usage()
{
  bigint bytes = AtomVecSphere::memory_usage();

  if (atom->memcheck("outerMass")) bytes += memory->usage(outerMass,nmax);
  if (atom->memcheck("outerRadius")) bytes += memory->usage(outerRadius,nmax);
 // if (atom->memcheck("atom_q")) bytes += memory->usage(atom_q,nmax);

  return bytes;
}

void AtomVecBio::set_group_mask() {

  for (int i = 1; i < group->ngroup; i++) {
    if (strcmp(group->names[i],"EPS") == 0) {
      maskEPS = pow(2, i) + 1;
    } else if (strcmp(group->names[i],"HET") == 0) {
      maskHET = pow(2, i) + 1;
    } else if (strcmp(group->names[i],"DEAD") == 0) {
      maskDEAD = pow(2, i) + 1;
    }
  }
}

