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

#include "compute_plasmid_ave_plasmid.h"

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
ComputePlasmidAvePlasmid::ComputePlasmidAvePlasmid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal nufeb/plasmid/nbirth command");
  fix_plasmid = nullptr;

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for compute plasmid/ave_plasmid");
  fix_plasmid = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  vector_flag = 1;
  extvector = 0;
  size_vector = 2;

  memory->create(vector,size_vector,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputePlasmidAvePlasmid::~ComputePlasmidAvePlasmid()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputePlasmidAvePlasmid::compute_vector()
{
  invoked_vector = update->ntimestep;
  int *mask = atom->mask;
  double ave_plm = 0.0;
  int nbacilli = 0;
  double sd = 0.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      ave_plm += static_cast<int>(fix_plasmid->vprop[i]);
      nbacilli++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &ave_plm, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &nbacilli, 1, MPI_INT, MPI_SUM, world);

  if (nbacilli) ave_plm = (double)ave_plm/nbacilli;

  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      sd += (static_cast<int>(fix_plasmid->vprop[i]) - ave_plm)*
	  (static_cast<int>(fix_plasmid->vprop[i]) - ave_plm);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &sd, 1, MPI_DOUBLE, MPI_SUM, world);

  if (nbacilli) sd /= nbacilli;

  vector[0] = ave_plm;
  vector[1] = sqrt(sd);
}
