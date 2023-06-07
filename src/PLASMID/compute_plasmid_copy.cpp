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

#include "compute_plasmid_copy.h"

#include "atom.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "domain.h"
#include "memory.h"
#include "math_const.h"
#include "atom_vec_bacillus.h"
#include "group.h"
#include <stdio.h>

#include "comm.h"
#include "fix_property_plasmid.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
ComputePlasmidCopy::ComputePlasmidCopy(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute nufeb/plasmid/copy command");
  avec = nullptr;
  fix_plasmid = nullptr;

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"compute nufeb/plasmid/copy requires "
      "atom style bacillus");

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for compute plasmid/copy");
  fix_plasmid = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  vector_flag = 1;
  size_vector = utils::inumeric(FLERR,arg[3],true,lmp);
  vector = new double[size_vector];

  if (narg != size_vector+4)
    error->all(FLERR,"Illegal compute nufeb/plasmid/copy command");

  cp = memory->create(cp, size_vector, "compute nufeb/plasmid/copy: cp");
  cpflag = memory->create(cpflag, size_vector, "compute nufeb/plasmid/copy: cpflag");

  for (int i = 0; i < size_vector; i++) {
    int n = strlen(arg[4+i]);
    char *str;

    if (arg[4+i][n-1] == '+') {
      str = new char[n];
      cpflag[i] = 1;
      strncpy(str, arg[4+i], n-1);
      str[n-1] = '\0';
    } else if (arg[4+i][n-1] == '-') {
      str = new char[n];
      cpflag[i] = 2;
      strncpy(str, arg[4+i], n-1);
      str[n-1] = '\0';
    } else if (isdigit(arg[4+i][n-1])) {
      str = new char[n+1];
      cpflag[i] = 0;
      strcpy(str, arg[4+i]);
    } else
      error->all(FLERR,"Illegal compute nufeb/plasmid/copy command");

    cp[i] = utils::inumeric(FLERR,str,true,lmp);
    vector[i] = 0;

    delete[] str;
  }

}

/* ---------------------------------------------------------------------- */

ComputePlasmidCopy::~ComputePlasmidCopy()
{
  delete [] vector;
  memory->destroy(cp);
  memory->destroy(cpflag);
}

/* ---------------------------------------------------------------------- */

void ComputePlasmidCopy::compute_vector()
{
  invoked_vector = update->ntimestep;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < size_vector; i++) {
    vector[i] = 0;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (int j = 0; j < size_vector; j++) {
	if (cpflag[j] == 0 && static_cast<int>(fix_plasmid->vprop[i]) == cp[j])
	  vector[j]++;
	if (cpflag[j] == 1 && static_cast<int>(fix_plasmid->vprop[i]) >= cp[j])
	  vector[j]++;
	if (cpflag[j] == 2 && static_cast<int>(fix_plasmid->vprop[i]) <= cp[j])
	  vector[j]++;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,vector,size_vector,MPI_DOUBLE,MPI_SUM,world);
}

