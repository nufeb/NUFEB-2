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

#include "compute_plasmid_nbirth.h"

#include <cstring>

#include "atom.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "domain.h"
#include "math_const.h"
#include "atom_vec_bacillus.h"
#include "fix_property_plasmid.h"
#include "fix_store.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
ComputePlasmidNBirth::ComputePlasmidNBirth(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal nufeb/plasmid/nbirth command");
  fix_plasmid = nullptr;

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for compute plasmid/nbirth");
  fix_plasmid = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  scalar = 0.0;
  nbacilli = 0;
  nbirth = 0;
  flag = 0;

  vector_flag = 1;
  extvector = 0;
  size_vector = 2;

  memory->create(vector,size_vector,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputePlasmidNBirth::~ComputePlasmidNBirth()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputePlasmidNBirth::compute_vector()
{
  invoked_vector = update->ntimestep;
  double ave_nbirth = 0.0;
  double sd = 0.0;

 // MPI_Allreduce(MPI_IN_PLACE, &nbirth, 1, MPI_DOUBLE, MPI_SUM, world);
 // MPI_Allreduce(MPI_IN_PLACE, &nbacilli, 1, MPI_INT, MPI_SUM, world);

  if (nbacilli) ave_nbirth = (double)nbirth/nbacilli;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (flag)sd += (static_cast<int>(fix_plasmid->vprop[i]) - ave_nbirth)*
	  (static_cast<int>(fix_plasmid->vprop[i]) - ave_nbirth);
    }
  }

  //MPI_Allreduce(MPI_IN_PLACE, &sd, 1, MPI_DOUBLE, MPI_SUM, world);

  if (nbacilli) sd /= nbacilli;

  vector[0] = ave_nbirth;
  vector[1] = sqrt(sd);
  flag = 0;
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputePlasmidNBirth::set_arrays(int i)
{
  flag = 1;
  if (atom->mask[i] & groupbit) {
    nbacilli++;
    nbirth += static_cast<int>(fix_plasmid->vprop[i]);
  }
}
