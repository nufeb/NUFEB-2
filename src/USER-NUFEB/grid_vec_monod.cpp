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

#include "grid_vec_monod.h"
#include "grid.h"
#include "error.h"
#include "memory.h"
#include "grid_masks.h"
#include "comm.h"
#include "atom.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecMonod::GridVecMonod(LAMMPS *lmp) : GridVec(lmp)
{
  mask = NULL;
  conc = NULL;
  reac = NULL;
  dens = NULL;
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  if (n > nmax) {
    mask = memory->grow(grid->mask, n, "nufeb/monod:mask");
    conc = memory->grow(grid->conc, grid->nsubs, n, "nufeb/monod:conc");
    reac = memory->grow(grid->reac, grid->nsubs, n, "nufeb/monod:reac");
    dens = memory->grow(grid->dens, group->ngroup, n, "nufeb/monod:dens");
    growth = memory->grow(grid->growth, group->ngroup, n, 2, "nufeb/monod:grow");
    nmax = n;
    grid->nmax = nmax;

    grid->mask = mask;
    grid->conc = conc;
    grid->reac = reac;
    grid->dens = dens;
    grid->growth = growth;
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonod::pack_comm(int n, int *cells, double *buf)
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

void GridVecMonod::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonod::pack_exchange(int n, int *cells, double *buf)
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

void GridVecMonod::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonod::set(int sub, double domain, double nx, double px,
		       double ny, double py, double nz, double pz)
{
  for (int i = 0; i < grid->ncells; i++) {
    if (!(mask[i] & CORNER_MASK)) {
      if (mask[i] & X_NB_MASK) {
	conc[sub][i] = nx;
      } else if (mask[i] & X_PB_MASK) {
	conc[sub][i] = px;
      } else if (mask[i] & Y_NB_MASK) {
	conc[sub][i] = ny;
      } else if (mask[i] & Y_PB_MASK) {
	conc[sub][i] = py;
      } else if (mask[i] & Z_NB_MASK) {
	conc[sub][i] = nz;
      } else if (mask[i] & Z_PB_MASK) {
	conc[sub][i] = pz;
      } else {
	conc[sub][i] = domain;
      }
    }
  }
}
