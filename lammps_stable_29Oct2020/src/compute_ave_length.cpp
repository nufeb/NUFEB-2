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

#include <math.h>

#include "compute_ave_length.h"
#include "atom.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "memory.h"
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

  vector_flag = 1;
  extvector = 0;
  size_vector = 2;

  memory->create(vector,size_vector,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputeAveLength::~ComputeAveLength()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputeAveLength::compute_vector()
{
  invoked_vector = update->ntimestep;
  int *mask = atom->mask;
  int nbacilli = 0;
  double sd = 0.0;
  double ave_len = 0.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (atom->bacillus_flag) {
	int ibonus = atom->bacillus[i];
	AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

	ave_len += bonus->length;
	nbacilli++;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &ave_len, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &nbacilli, 1, MPI_INT, MPI_SUM, world);

  if (nbacilli) ave_len /= nbacilli;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (atom->bacillus_flag) {
	int ibonus = atom->bacillus[i];
	AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

	sd += (bonus->length - ave_len)*(bonus->length - ave_len);
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &sd, 1, MPI_DOUBLE, MPI_SUM, world);

  if (nbacilli) sd /= nbacilli;

  vector[0] = ave_len;
  vector[1] = sqrt(sd);
}
