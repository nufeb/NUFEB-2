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

#include "grid_kokkos.h"
#include "grid_vec.h"
#include "grid_vec_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridKokkos::GridKokkos(LAMMPS *lmp) : Grid(lmp)
{
  gridKK = this;
}

/* ---------------------------------------------------------------------- */

GridKokkos::~GridKokkos()
{
  memoryKK->destroy_kokkos(k_mask, mask);
  memoryKK->destroy_kokkos(k_bulk, bulk);
  memoryKK->destroy_kokkos(k_conc, conc);
  memoryKK->destroy_kokkos(k_reac, reac);
  memoryKK->destroy_kokkos(k_diff_coeff, diff_coeff);
  memoryKK->destroy_kokkos(k_dens, dens);
  memoryKK->destroy_kokkos(k_boundary, boundary);
  memoryKK->destroy_kokkos(k_growth, growth);
}

/* ---------------------------------------------------------------------- */

void GridKokkos::sync(const ExecutionSpace space, unsigned int mask)
{
  if (space == Device && lmp->kokkos->auto_sync)
    ((GridVecKokkos *) gvec)->modified(Host,mask);

  ((GridVecKokkos *) gvec)->sync(space,mask);
}

/* ---------------------------------------------------------------------- */

void GridKokkos::modified(const ExecutionSpace space, unsigned int mask)
{
  ((GridVecKokkos *) gvec)->modified(space,mask);

  if (space == Device && lmp->kokkos->auto_sync)
    ((GridVecKokkos *) gvec)->sync(Host,mask);
}

/* ---------------------------------------------------------------------- */

void GridKokkos::sync_overlapping_device(const ExecutionSpace space, unsigned int mask)
{
  ((GridVecKokkos *) gvec)->sync_overlapping_device(space,mask);
}

/* ---------------------------------------------------------------------- */

void GridKokkos::sync_modify(ExecutionSpace execution_space,
                             unsigned int datamask_read,
                             unsigned int datamask_modify)
{
  sync(execution_space,datamask_read);
  modified(execution_space,datamask_modify);
}

/* ---------------------------------------------------------------------- */

GridVec *GridKokkos::new_gvec(const char *style, int trysuffix, int &sflag)
{
  GridVec* gvec = Grid::new_gvec(style,trysuffix,sflag);
  if (!gvec->kokkosable)
    error->all(FLERR,"KOKKOS package requires a kokkos enabled grid_style");
  return gvec;
}
