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

#include <string.h>
#include <math.h>
#include "atom.h"
#include "error.h"
#include "modify.h"
#include "group.h"
#include "update.h"
#include "memory.h"
#include "grid.h"
#include "atom_masks.h"
#include "fix_plasmid_conjugation.h"
#include "fix_property_plasmid.h"

#include "fix.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPlasmidConjugation::FixPlasmidConjugation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), pilus(nullptr), tpilus(nullptr)
{
  if (narg < 5)
    error->all(FLERR, "Illegal nufeb/plasmid/replication command");

  fix_plm = nullptr;
  k_pilus = 1e-5;
  npili = 0;
  pilus_len = 8e-9;

  seed = utils::inumeric(FLERR,arg[3],true,lmp);
  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  int ifix = modify->find_fix_by_style("^nufeb/property/plasmid");
  if (ifix < 0 ) error->all(FLERR,"Illegal nufeb/plasmid/replication command: "
      "requires fix nufeb/property/plasmid");
  fix_plm = (FixPropertyPlasmid *) modify->fix[ifix];

  irecipient = group->find(arg[4]);
  itrans = group->find(arg[5]);

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "length") == 0) {
      pilus_len = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "init") == 0) {
      plm_init = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      if (plm_init < 0)
	error->all(FLERR,"Illegal fix nufeb/property/plasmid command: init");
      iarg += 2;
    }  else {
      error->all(FLERR,"Illegal fix nufeb/plasmid/replication command");
    }
  }

  fix_plm->rep_flag = 1;
  grow_arrays(atom->nmax);
}


/* ---------------------------------------------------------------------- */
FixPlasmidConjugation:: ~FixPlasmidConjugation() {
  memory->destroy(pilus);
  memory->destroy(tpilus);

  delete random;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPlasmidConjugation::grow_arrays(int nmax)
{
  pili_max = nmax / 2;
  memory->grow(pilus,pili_max,2,"fix_nufeb/plasmid/conjugate:pilus");
  memory->grow(tpilus,pili_max,"fix_nufeb/plasmid/conjugate:pili_max");
}


/* ---------------------------------------------------------------------- */

int FixPlasmidConjugation::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPlasmidConjugation::biology_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixPlasmidConjugation::compute()
{

}
