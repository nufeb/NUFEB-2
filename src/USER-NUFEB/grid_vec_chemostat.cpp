/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "grid_vec_chemostat.h"

#include "grid.h"
#include "force.h"
#include "error.h"
#include "memory.h"
#include "grid_masks.h"
#include "comm.h"
#include "atom.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecChemostat::GridVecChemostat(LAMMPS *lmp) : GridVec(lmp)
{
  mask = nullptr;
  conc = nullptr;
  reac = nullptr;
  dens = nullptr;
  growth = nullptr;
  bulk = nullptr;
  mw = nullptr;
  boundary = nullptr;
  diff_coeff = nullptr;
  grid->chemostat_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GridVecChemostat::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecChemostat::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  if (n > nmax) {
    mask = memory->grow(grid->mask, n, "nufeb/monod:mask");
    conc = memory->grow(grid->conc, grid->nsubs, n, "nufeb/chemostat:conc");
    reac = memory->grow(grid->reac, grid->nsubs, n, "nufeb/chemostat:reac");
    dens = memory->grow(grid->dens, group->ngroup, n, "nufeb/chemostat:dens");
    growth = memory->grow(grid->growth, group->ngroup, n, 2, "nufeb/chemostat:grow");
    boundary = memory->grow(grid->boundary, grid->nsubs, 6, "nufeb/chemostat:boundary");
    diff_coeff = memory->grow(grid->diff_coeff, grid->nsubs, n, "nufeb/chemostat:diff_coeff");
    bulk = memory->grow(grid->bulk, grid->nsubs, "nufeb/chemostat:bulk");
    mw = memory->grow(grid->mw, grid->nsubs, "nufeb/chemostat:mw");

    nmax = n;
    grid->nmax = nmax;

    grid->mask = mask;
    grid->conc = conc;
    grid->reac = reac;
    grid->dens = dens;
    grid->growth = growth;
    grid->bulk = bulk;
    grid->mw = mw;
    grid->boundary = boundary;
    grid->diff_coeff = diff_coeff;
  }
}

/* ---------------------------------------------------------------------- */

int GridVecChemostat::pack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecChemostat::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int GridVecChemostat::pack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecChemostat::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridVecChemostat::set(int narg, char **arg)
{
  if (narg < 6) error->all(FLERR, "Invalid grid_modify set command");
  int isub = grid->find(arg[1]);
  if (isub < 0) error->all(FLERR,"Cannot find substrate name");

  for (int i = 0; i < 3; i++) {
    if ((arg[2+i][0] == 'p' && arg[2+i][1] != 'p') ||
	(arg[2+i][1] == 'p' && arg[2+i][0] != 'p'))
      error->all(FLERR, "Illegal boundary condition: unpaired periodic BC");

    for (int j = 0; j < 2; j++) {
      if (arg[2+i][j] == 'p') {
        boundary[isub][2*i+j] = PERIODIC;
        grid->periodic[i] = 1;
      } else if (arg[2+i][j] == 'n') {
        boundary[isub][2*i+j] = NEUMANN;
      } else if (arg[2+i][j] == 'd') {
        boundary[isub][2*i+j] = DIRICHLET;
      } else {
        error->all(FLERR, "Illegal boundary condition: unknown keyword");
      }
    }
  }

  double domain = utils::numeric(FLERR,arg[5],true,lmp);
  if (domain < 0) error->all(FLERR, "Illegal initial substrate concentration");

  bulk[isub] = domain;
  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "bulk") == 0) {
      bulk[isub] = utils::numeric(FLERR, arg[iarg+1], true, lmp);
      if (bulk[isub] < 0) error->all(FLERR, "Illegal initial bulk concentration");
      iarg += 2;
    } else if (strcmp(arg[iarg], "mw") == 0) {
      mw[isub] = utils::numeric(FLERR, arg[iarg+1], true, lmp);
      if (mw[isub] <= 0) error->all(FLERR, "Illegal molecular weight value");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal grid_modify command");
    }
  }

  set_grid(isub, domain, bulk[isub]);
}


/* ---------------------------------------------------------------------- */

void GridVecChemostat::set_grid(int isub, double domain, double bulk)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK)) {
      if (((mask[i] & X_NB_MASK) && (boundary[isub][0] == DIRICHLET)) ||
	  ((mask[i] & X_PB_MASK) && (boundary[isub][1] == DIRICHLET)) ||
	  ((mask[i] & Y_NB_MASK) && (boundary[isub][2] == DIRICHLET)) ||
	  ((mask[i] & Y_PB_MASK) && (boundary[isub][3] == DIRICHLET)) ||
	  ((mask[i] & Z_NB_MASK) && (boundary[isub][4] == DIRICHLET)) ||
	  ((mask[i] & Z_PB_MASK) && (boundary[isub][5] == DIRICHLET))) {
        conc[isub][i] = bulk;
      } else {
    	conc[isub][i] = domain;
      }
    }
    grid->reac[isub][i] = 0.0;
    grid->diff_coeff[isub][i] = 0.0;
  }
}
