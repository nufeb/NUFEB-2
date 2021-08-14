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

#include "compute_ave_length.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "math_const.h"
#include "atom_vec_bacillus.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeAveLength::ComputeAveLength(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute nufeb/volume command");
  avec = nullptr;
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"compute nufeb/plasmid requires "
      "atom style bacillus");

  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputeAveLength::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  int *mask = atom->mask;

  scalar = 0.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((mask[i] & groupbit)) {
      if (atom->bacillus_flag) {
	int ibonus = atom->bacillus[i];
	AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

	scalar += bonus->length;
      }
    }
  }
  if (atom->nlocal)
    scalar /= atom->nlocal;

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world); 

  return scalar/comm->nprocs;
}
