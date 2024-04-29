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

#include "fix_property_ancestor.h"

#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAncestor::FixPropertyAncestor(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix property/nufeb/ancestor command");

  create_attribute = 1;
  // use vprop if size_peratom_cols = 1
  size_peratom_cols = 1;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

// The first atoms have are their own ancestor, so 
// set their vprop to their identity (tag)
void FixPropertyAncestor::init()
{
  //printf("Initializing atoms:\n"); 
  for (int i = 0; i < atom->nlocal; i++) {
      ///printf("\tAtom i: %d tag: %d\n",i,atom->tag[i]);
      vprop[i] = atom->tag[i];
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyAncestor::set_arrays(int j)
{
    vprop[j] = -1; // no ancestor assigned
}

/* ----------------------------------------------------------------------
   update array values of two daughter cells i, j
   called in fix_divide
------------------------------------------------------------------------- */
void FixPropertyAncestor::update_arrays(int i, int j)
{
  // The newly created cell should have the same ancestor as its sister
  vprop[j] = vprop[i];
  //copy_arrays(i,j,0);  //didn't make a difference
}

int FixPropertyAncestor::setmask()
{
  return 0;
}


