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

#include "compute_bio_diversity.h"

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebDiversity::ComputeNufebDiversity(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ntypes command");

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebDiversity::~ComputeNufebDiversity()
{
}

/* ---------------------------------------------------------------------- */

double ComputeNufebDiversity::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  int *ntypeList = new int[ntypes + 1]();

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      ntypeList[t] ++;
    }
  }

  int n = 0;
  for (int i = 1; i < ntypes + 1; i++) {
    n = n + ntypeList[i] * (ntypeList[i] - 1);
  }

  int m = nlocal * (nlocal - 1);
  double x = n/(double)m;
  scalar = 1 - x;

  delete [] ntypeList;

  return scalar;
}
