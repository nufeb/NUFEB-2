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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebDiversity::ComputeNufebDiversity(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ntypes command");

  vector_flag = 1;
  size_vector = atom->ntypes;
  extvector = 0;

  vector = new double[atom->ntypes]();
}

/* ---------------------------------------------------------------------- */

ComputeNufebDiversity::~ComputeNufebDiversity()
{
   delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeNufebDiversity::compute_vector()
{
  invoked_vector = update->ntimestep;

  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  int *ntypeList = new int[ntypes + 1]();

  for (int i = 0; i < ntypes; i++) {
    vector[i] = 0;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      int t = type[i];
      ntypeList[t] ++;
    }
  }

  for (int i = 1; i < ntypes + 1; i++) {
    int sum = 0;
    for (int j = 1; j <= ntypeList[i]; j++) {
      sum = sum + j * (j - 1);
    }

    int a = nlocal * (nlocal - 1);

    if (a != 0)
      vector[i-1] = 1 - ((double)sum / a);
    else
      vector[i-1] = 1;
  }

  delete [] ntypeList;
}
