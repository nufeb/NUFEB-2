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
  avec = nullptr;
  fix_plasmid = nullptr;;

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"compute nufeb/plasmid/nbirth requires "
      "atom style bacillus");

  int ifix = modify->find_fix_by_style("^nufeb/property/plasmid");
  if (ifix < 0 ) error->all(FLERR,"Illegal nufeb/plasmid/nbirth command: requires fix nufeb/property/plasmid");
  fix_plasmid = (FixPropertyPlasmid *) modify->fix[ifix];

  scalar = 0.0;
  ncell = 0;
  nplm = 0;

  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

double ComputePlasmidNBirth::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  if (ncell)
    scalar = (double)nplm/ncell;

  return scalar;
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputePlasmidNBirth::set_arrays(int i)
{
  ncell++;
  nplm += (int)fix_plasmid->vprop[i];
}
