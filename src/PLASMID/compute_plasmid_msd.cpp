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

#include "compute_plasmid_msd.h"

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
ComputePlasmidMSD::ComputePlasmidMSD(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_fix(nullptr)
{
  if (narg < 3) error->all(FLERR,"Illegal compute nufeb/plasmid/msd command");
  avec = nullptr;
  fix_plasmid = nullptr;

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"compute nufeb/plasmid/msd requires "
      "atom style bacillus");

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for compute/plasmid/msd");
  fix_plasmid = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  vector_flag = 1;
  size_vector = 4;

  nmsd = 0;

  // create a new fix STORE style for reference positions
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int ncol = fix_plasmid->plm_max*3;
  std::string fixcmd = id + std::string("_COMPUTE_STORE");
  id_fix = new char[fixcmd.size()+1];
  strcpy(id_fix,fixcmd.c_str());

  fixcmd += fmt::format(" {} STORE peratom 1 {}",group->names[igroup], ncol);

  modify->add_fix(fixcmd);
  fix_store = (FixStore *) modify->fix[modify->nfix-1];
  // displacement vector

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

void ComputePlasmidMSD::init()
{
  vector[0] = vector[1] = vector[2] = vector[3] = 0.0;
  // calculate x,y,z for fix store array
  // skip if reset from restart file

  if (fix_store->restart_reset) fix_store->restart_reset = 0;
  else {
    double **xpm_original = fix_store->astore;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	for (int j = 0; j < (int)fix_plasmid->vprop[i]; j++) {
	  int jx0 = j*3;
	  int jx1 = j*3+1;
	  int jx2 = j*3+2;

	  xpm_original[i][jx0] = fix_plasmid->plm_x[i][jx0];
	  xpm_original[i][jx1] = fix_plasmid->plm_x[i][jx1];
	  xpm_original[i][jx2] = fix_plasmid->plm_x[i][jx2];
	}
      }
    }
  }

  // set fix which stores reference atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute msd fix ID");
  fix_store = (FixStore *) modify->fix[ifix];

}

/* ---------------------------------------------------------------------- */

ComputePlasmidMSD::~ComputePlasmidMSD()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputePlasmidMSD::compute_vector()
{

  invoked_vector = update->ntimestep;

  double **xpm_original = fix_store->astore;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dx,dy,dz;

  double msd[4];
  msd[0] = msd[1] = msd[2] = msd[3] = 0.0;
  nmsd = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      for (int j = 0; j < (int)fix_plasmid->vprop[i]; j++) {
	int jx0 = j*3;
	int jx1 = j*3+1;
	int jx2 = j*3+2;

	dx = fix_plasmid->plm_x[i][jx0] - xpm_original[i][jx0];
	dy = fix_plasmid->plm_x[i][jx1] - xpm_original[i][jx1];
	dz = fix_plasmid->plm_x[i][jx2] - xpm_original[i][jx2];

	msd[0] += dx*dx;
	msd[1] += dy*dy;
	msd[2] += dz*dz;
	msd[3] += dx*dx + dy*dy + dz*dz;
	// nmsd = # of plasmid in group
	nmsd++;
      }
    }

  MPI_Allreduce(msd,vector,4,MPI_DOUBLE,MPI_SUM,world);

  if (nmsd) {
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] /= nmsd;
    vector[3] /= nmsd;
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's storage values, called when atom is created
------------------------------------------------------------------------- */

void ComputePlasmidMSD::set_arrays(int i)
{
  double **xpm_original = fix_store->astore;

  for (int j = 0; j < (int)fix_plasmid->vprop[i]; j++) {
    int jx0 = j*3;
    int jx1 = j*3+1;
    int jx2 = j*3+2;

    xpm_original[i][jx0] = fix_plasmid->plm_x[i][jx0];
    xpm_original[i][jx1] = fix_plasmid->plm_x[i][jx1];
    xpm_original[i][jx2] = fix_plasmid->plm_x[i][jx2];
  }
}
