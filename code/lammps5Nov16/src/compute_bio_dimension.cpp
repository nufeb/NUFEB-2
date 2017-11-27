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

#include "compute_bio_dimension.h"

#include <math.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNufebDimension::ComputeNufebDimension(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute fractal dimension command");

  scalar_flag = 1;
  extscalar = 0;
}

/* ---------------------------------------------------------------------- */

ComputeNufebDimension::~ComputeNufebDimension()
{
}

/* ---------------------------------------------------------------------- */
double ComputeNufebDimension::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  scalar = 0;

  double md = 0.0;
  double m = 0.0;
  double r = 0.0;

  double ra = 0.0;
  double rm = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double dia = radius[i] * 2;
      md += rmass[i] * dia * dia;
      m += rmass[i];
      r += radius[i];
    }
  }

  ra = pow(md/m, 1.0/2);
  rm = r / nlocal;

  scalar = log(ra/m) / log(nlocal);

  return scalar;
}
