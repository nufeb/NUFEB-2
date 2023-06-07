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

#include "grid_vec_simple.h"

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

GridVecSimple::GridVecSimple(LAMMPS *lmp) : GridVec(lmp)
{
  mask = nullptr;
  conc = nullptr;

  grid->simple_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GridVecSimple::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecSimple::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  if (n > nmax) {
    mask = memory->grow(grid->mask, n, "nufeb/monod:mask");
    conc = memory->grow(grid->conc, grid->nsubs, n, "nufeb/chemostat:conc");

    nmax = n;
    grid->nmax = nmax;

    grid->mask = mask;
    grid->conc = conc;
  }
}

/* ---------------------------------------------------------------------- */

int GridVecSimple::pack_comm(int n, int *cells, double *buf)
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

void GridVecSimple::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int GridVecSimple::pack_exchange(int n, int *cells, double *buf)
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

void GridVecSimple::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridVecSimple::set(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR, "Invalid grid_modify set command");
  int isub = grid->find(arg[1]);
  if (isub < 0) error->all(FLERR,"Cannot find substrate name");

  double domain = utils::numeric(FLERR,arg[2],true,lmp);
  if (domain < 0) error->all(FLERR, "Illegal initial substrate concentration");

  set_grid(isub, domain);
}


/* ---------------------------------------------------------------------- */

void GridVecSimple::set_grid(int isub, double domain)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK)) {
      conc[isub][i] = domain;
    }
  }
}
