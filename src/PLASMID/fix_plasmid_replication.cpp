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
#include "fix_plasmid_replication.h"
#include "fix_property_plasmid.h"

#include "fix.h"
#include "random_park.h"
#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define DELTA 1.005
/* ---------------------------------------------------------------------- */

FixPlasmidReplication::FixPlasmidReplication(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR, "Illegal nufeb/plasmid/replication command");

  fix_plm = nullptr;
  nproteins = nullptr;

  mean_protein = 20;
  init_protein = 0;
  alpha = 1.0;

  seed = utils::inumeric(FLERR,arg[3],true,lmp);
  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for fix plasmid/replicatio");
  fix_plm = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  int iarg = 4;
  while (iarg < narg) {
     if (strcmp(arg[iarg], "mean") == 0) {
      mean_protein = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "init") == 0) {
      init_protein = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "alpha") == 0) {
      alpha = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/plasmid/replication command");
    }
  }

  fix_plm->rep_flag = 1;
}


/* ---------------------------------------------------------------------- */
FixPlasmidReplication:: ~FixPlasmidReplication() {
  memory->destroy(nproteins);
  delete random;
}

/* ---------------------------------------------------------------------- */
void FixPlasmidReplication::init() {
  grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nlocal; i++) {
    // initialise protein number
    for (int j = 0; j < fix_plm->plm_max; j++) {
      nproteins[i][j] = init_protein;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPlasmidReplication::grow_arrays(int nmax)
{
  memory->grow(nproteins,nmax,fix_plm->plm_max,"fix_nufeb/plasmid/replication:nproteins");
  fix_plm->nproteins = nproteins;
}


/* ---------------------------------------------------------------------- */

int FixPlasmidReplication::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPlasmidReplication::biology_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixPlasmidReplication::compute()
{
  for (int i = 0; i < atom->nlocal; i++){
    if (atom->mask[i] & groupbit) {
      replication(i);
    }
  }
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */
void FixPlasmidReplication::replication(int i)
{
  double current_t = update->ntimestep * update->dt;
  double next_t = (update->ntimestep + 1) * update->dt;

  const int cell = grid->cell(atom->x[i]);
  double growth = grid->growth[igroup][cell][0];

  double *plm_x = fix_plm->plm_x[i];
  double *vprop = fix_plm->vprop;
  double plm_dia = fix_plm->plm_dia;

  // update initiator proteins
  for (int m = 0; m < static_cast<int>(vprop[i]); m++)
    nproteins[i][m] += alpha * growth * atom->rmass[i] * update->dt;

  for (int m = 0; m < static_cast<int>(vprop[i]); m++) {
    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;
    double max = mean_protein + mean_protein*0.1*random->gaussian();

    if (nproteins[i][m] > max && vprop[i] < fix_plm->plm_max) {
      double ilimit[3];

      int n = static_cast<int>(vprop[i]);
      int n0 = n*3;
      int n1 = n*3+1;
      int n2 = n*3+2;

      double theta = random->uniform() * 2 * MY_PI;
      double phi = random->uniform() * (MY_PI);

      plm_x[n0] = plm_x[m0] + (plm_dia * cos(theta) * sin(phi) * DELTA);
      plm_x[n1] = plm_x[m1] + (plm_dia * sin(theta) * sin(phi) * DELTA);
      plm_x[n2] = plm_x[m2] + (plm_dia * cos(phi) * DELTA);

      fix_plm->relocate_plm_x(i,m);
      vprop[i]++;

      nproteins[i][m] = 0.0;
      nproteins[i][n] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPlasmidReplication::memory_usage()
{
  double bytes;

  bytes += atom->nmax*fix_plm->plm_max*sizeof(double);

  return bytes;
}
