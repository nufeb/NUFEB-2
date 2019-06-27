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

#include "grid_vec_monod_kokkos.h"
#include "grid_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "grid_masks.h"
#include "comm.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecMonodKokkos::GridVecMonodKokkos(LAMMPS *lmp) : GridVecKokkos(lmp)
{
  mask = NULL;
  conc = NULL;
  reac = NULL;
  dens = NULL;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");
  
  if (n > nmax) {
    sync(Device, ALL_MASK);
    modified(Device, ALL_MASK);

    memoryKK->grow_kokkos(gridKK->k_mask, gridKK->mask, n, "nufeb/monod:mask");
    memoryKK->grow_kokkos(gridKK->k_conc, gridKK->conc, grid->nsubs, n, "nufeb/monod:conc");
    memoryKK->grow_kokkos(gridKK->k_reac, gridKK->reac, grid->nsubs, n, "nufeb/monod:reac");
    memoryKK->grow_kokkos(gridKK->k_dens, gridKK->dens, group->ngroup, n, "nufeb/monod:dens");
    memoryKK->grow_kokkos(gridKK->k_growth, gridKK->growth, group->ngroup, n, 2, "nufeb/monod:grow");
    nmax = n;
    grid->nmax = nmax;

    mask = gridKK->mask;
    d_mask = gridKK->k_mask.d_view;
    h_mask = gridKK->k_mask.h_view;
    conc = gridKK->conc;
    d_conc = gridKK->k_conc.d_view;
    h_conc = gridKK->k_conc.h_view;
    reac = gridKK->reac;
    d_reac = gridKK->k_reac.d_view;
    h_reac = gridKK->k_reac.h_view;
    dens = gridKK->dens;
    d_dens = gridKK->k_dens.d_view;
    h_dens = gridKK->k_dens.h_view;
    growth = gridKK->growth;
    d_growth = gridKK->k_growth.d_view;
    h_growth = gridKK->k_growth.h_view;

    sync(Host, ALL_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_comm(int n, int *cells, double *buf)
{
  sync(Host, CONC_MASK);
  
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_comm(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }

  modified(Host, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_exchange(int n, int *cells, double *buf)
{
  sync(Host, CONC_MASK);
  
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      buf[m++] = conc[s][cells[c]];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_exchange(int n, int *cells, double *buf)
{
  int m = 0;
  for (int s = 0; s < grid->nsubs; s++) {
    for (int c = 0; c < n; c++) {
      conc[s][cells[c]] = buf[m++];
    }
  }

  modified(Host, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_comm_kokkos(int, const DAT::tdual_int_2d &list, const DAT::tdual_xfloat_2d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_comm_kokkos(int, const DAT::tdual_int_2d &list, const DAT::tdual_xfloat_2d &buf)
{

}

/* ---------------------------------------------------------------------- */

int GridVecMonodKokkos::pack_exchange_kokkos(int, const DAT::tdual_int_2d &list, const DAT::tdual_xfloat_2d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::unpack_exchange_kokkos(int, const DAT::tdual_int_2d &list, const DAT::tdual_xfloat_2d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::set(int sub, double domain, double nx, double px,
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

  modified(Host, CONC_MASK);
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & GMASK_MASK) gridKK->k_mask.sync<LMPDeviceType>();
    if (mask & CONC_MASK) gridKK->k_conc.sync<LMPDeviceType>();
    if (mask & REAC_MASK) gridKK->k_reac.sync<LMPDeviceType>();
    if (mask & DENS_MASK) gridKK->k_dens.sync<LMPDeviceType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.sync<LMPDeviceType>();
  } else {
    if (mask & GMASK_MASK) gridKK->k_mask.sync<LMPHostType>();
    if (mask & CONC_MASK) gridKK->k_conc.sync<LMPHostType>();
    if (mask & REAC_MASK) gridKK->k_reac.sync<LMPHostType>();
    if (mask & DENS_MASK) gridKK->k_dens.sync<LMPHostType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & GMASK_MASK) gridKK->k_mask.modify<LMPDeviceType>();
    if (mask & CONC_MASK) gridKK->k_conc.modify<LMPDeviceType>();
    if (mask & REAC_MASK) gridKK->k_reac.modify<LMPDeviceType>();
    if (mask & DENS_MASK) gridKK->k_dens.modify<LMPDeviceType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.modify<LMPDeviceType>();
  } else {
    if (mask & GMASK_MASK) gridKK->k_mask.modify<LMPHostType>();
    if (mask & CONC_MASK) gridKK->k_conc.modify<LMPHostType>();
    if (mask & REAC_MASK) gridKK->k_reac.modify<LMPHostType>();
    if (mask & DENS_MASK) gridKK->k_dens.modify<LMPHostType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void GridVecMonodKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if ((mask & GMASK_MASK) && gridKK->k_mask.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_1d>(gridKK->k_mask, space);
    if ((mask & CONC_MASK) && gridKK->k_conc.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_conc, space);
    if ((mask & REAC_MASK) && gridKK->k_reac.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_reac, space);
    if ((mask & DENS_MASK) && gridKK->k_dens.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_dens, space);
    if ((mask & GROWTH_MASK) && gridKK->k_growth.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_3d>(gridKK->k_growth, space);
  } else {
    if ((mask & GMASK_MASK) && gridKK->k_mask.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_1d>(gridKK->k_mask, space);
    if ((mask & CONC_MASK) && gridKK->k_conc.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_conc, space);
    if ((mask & REAC_MASK) && gridKK->k_reac.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_reac, space);
    if ((mask & DENS_MASK) && gridKK->k_dens.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_dens, space);
    if ((mask & GROWTH_MASK) && gridKK->k_growth.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_3d>(gridKK->k_growth, space);
  }
}
