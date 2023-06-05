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

#include "fix_merge_atom.h"

#include <math.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "pair.h"
#include "update.h"
#include "memory.h"
#include "force.h"
#include "math_const.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixMergeAtom::FixMergeAtom(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg < 5)
    error->all(FLERR, "Illegal fix nufeb/merge_eps command");

  list = nullptr;
  eps_den = 30;

  max_dia = utils::numeric(FLERR,arg[3],true,lmp);
  if (max_dia < 0)
    error->all(FLERR, "Illegal parameter in fix nufeb/merge_eps");
  seed = utils::numeric(FLERR,arg[4],true,lmp);

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "epsdens") == 0) {
      eps_den = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    }
  }

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  force_reneighbor = 1;
}

/* ---------------------------------------------------------------------- */

FixMergeAtom::~FixMergeAtom()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixMergeAtom::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMergeAtom::init_list(int id, NeighList *ptr)
{
  this->list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixMergeAtom::init()
{
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix  = 1;
  neighbor->requests[irequest]->size = 1;
}

/* ---------------------------------------------------------------------- */

void FixMergeAtom::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
  compute();
  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */
void FixMergeAtom::post_neighbor()
{
  // reset reneighbour flag
  next_reneighbor = 0;
}

/* ---------------------------------------------------------------------- */

void FixMergeAtom::compute()
{
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *mask = atom->mask;
  double **x = atom->x;
  double *radius = atom->radius;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  AtomVec *avec = atom->avec;

  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;
  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));

  double delx,dely,delz,rsq,min;

  int *dlist;
  memory->create(dlist,nall,"merge_eps:dlist");
  for (int i = 0; i < nall; i++) dlist[i] = 0;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    if (!(mask[i] & groupbit) || dlist[i] == 1) continue;

    if (atom->radius[i] * 2 < max_dia) {
      // pick the closest neighbor to merge
      int pick = -1;
      int *jlist = firstneigh[i];
      int jnum = numneigh[i];
      min = MAXTAGINT;

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        if (!(mask[j] & groupbit) || dlist[j] == 1 || j >= nlocal) continue;

        delx = x[i][0]-x[j][0];
        dely = x[i][1]-x[j][1];
        delz = x[i][2]-x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < min) {
          min = rsq;
          pick = jj;
        }
      }

      if (pick < 0) continue;
      int j = jlist[pick];

      double density = atom->rmass[i] / (four_thirds_pi * atom->radius[i] * atom->radius[i] * atom->radius[i]);
      double new_mass = atom->rmass[i] + atom->rmass[j];
      double new_outer_mass = atom->outer_mass[i] + atom->outer_mass[j];
      double new_biomass = atom->biomass[i] + atom->biomass[j];
      double new_rad = pow(three_quarters_pi * new_mass / density, third);
      double new_outer_rad  = pow(three_quarters_pi * ((new_mass / density) + (new_outer_mass / eps_den)), third);

      // randomly choose atom i or j as new atom's attributes
      if (random->uniform() > 0.5) {
        atom->x[i][0] = atom->x[j][0];
        atom->x[i][1] = atom->x[j][1];
        atom->x[i][2] = atom->x[j][2];
        atom->v[i][0] = atom->v[j][0];
        atom->v[i][1] = atom->v[j][1];
        atom->v[i][2] = atom->v[j][2];
        atom->f[i][0] = atom->f[j][0];
        atom->f[i][1] = atom->f[j][1];
        atom->f[i][2] = atom->f[j][2];
        atom->omega[i][0] = atom->omega[j][0];
        atom->omega[i][1] = atom->omega[j][1];
        atom->omega[i][2] = atom->omega[j][2];
      }

      dlist[j] = 1;
      atom->rmass[i] = new_mass;
      atom->radius[i] = new_rad;
      atom->biomass[i] = new_biomass;
      atom->outer_mass[i] = new_outer_mass;
      atom->outer_radius[i] = new_outer_rad;
    }
  }

  int i = 0;
  while (i < nlocal) {
    if (dlist[i]) {
      avec->copy(nlocal-1,i,1);
      dlist[i] = dlist[nlocal-1];
      nlocal--;
    } else i++;
  }

  atom->nlocal = nlocal;
  memory->destroy(dlist);

  // reset atom->natoms
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped
  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}
