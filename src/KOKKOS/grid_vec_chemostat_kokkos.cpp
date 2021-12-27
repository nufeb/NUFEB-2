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

#include "grid_vec_chemostat_kokkos.h"
#include "grid_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "grid_masks.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "group.h"
#include "lammps.h"
#include "kokkos.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecChemostatKokkos::GridVecChemostatKokkos(LAMMPS *lmp) : GridVecKokkos(lmp)
{
  mask = NULL;
  bulk = NULL;
  conc = NULL;
  reac = NULL;
  diff_coeff = NULL;
  dens = NULL;
  boundary = NULL;
  growth = NULL;
  grid->chemostat_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::init()
{
  GridVec::init();

  size_forward = grid->nsubs;
  size_exchange = grid->nsubs;
}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::grow(int n)
{
  if (n < 0 || n > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");
  
  if (n > nmax) {
    sync(Device, ALL_MASK);
    modified(Device, ALL_MASK);

    memoryKK->grow_kokkos(gridKK->k_mask, gridKK->mask, n, "nufeb/chemostat:mask");
    memoryKK->grow_kokkos(gridKK->k_bulk, gridKK->bulk, grid->nsubs, "nufeb/chemostat:bulk");
    memoryKK->grow_kokkos(gridKK->k_conc, gridKK->conc, grid->nsubs, n, "nufeb/chemostat:conc");
    memoryKK->grow_kokkos(gridKK->k_reac, gridKK->reac, grid->nsubs, n, "nufeb/chemostat:reac");
    memoryKK->grow_kokkos(gridKK->k_diff_coeff, gridKK->diff_coeff, grid->nsubs, n, "nufeb/chemostat:diff_coeff");
    memoryKK->grow_kokkos(gridKK->k_dens, gridKK->dens, group->ngroup, n, "nufeb/chemostat:dens");
    memoryKK->grow_kokkos(gridKK->k_boundary, gridKK->boundary, grid->nsubs, 6, "nufeb/chemostat:boundary");
    memoryKK->grow_kokkos(gridKK->k_growth, gridKK->growth, group->ngroup, n, 2, "nufeb/chemostat:grow");
    nmax = n;
    grid->nmax = nmax;

    mask = gridKK->mask;
    d_mask = gridKK->k_mask.d_view;
    h_mask = gridKK->k_mask.h_view;
    bulk = gridKK->bulk;
    d_bulk = gridKK->k_bulk.d_view;
    h_bulk = gridKK->k_bulk.h_view;
    conc = gridKK->conc;
    d_conc = gridKK->k_conc.d_view;
    h_conc = gridKK->k_conc.h_view;
    reac = gridKK->reac;
    d_diff_coeff = gridKK->k_diff_coeff.d_view;
    h_diff_coeff = gridKK->k_diff_coeff.h_view;
    diff_coeff = gridKK->diff_coeff;
    d_reac = gridKK->k_reac.d_view;
    h_reac = gridKK->k_reac.h_view;
    dens = gridKK->dens;
    d_dens = gridKK->k_dens.d_view;
    h_dens = gridKK->k_dens.h_view;
    boundary = gridKK->boundary;
    d_boundary = gridKK->k_boundary.d_view;
    h_boundary = gridKK->k_boundary.h_view;
    growth = gridKK->growth;
    d_growth = gridKK->k_growth.d_view;
    h_growth = gridKK->k_growth.h_view;

    sync(Host, ALL_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int GridVecChemostatKokkos::pack_comm(int n, int *cells, double *buf)
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

void GridVecChemostatKokkos::unpack_comm(int n, int *cells, double *buf)
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

int GridVecChemostatKokkos::pack_exchange(int n, int *cells, double *buf)
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

void GridVecChemostatKokkos::unpack_exchange(int n, int *cells, double *buf)
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

template <class DeviceType>
struct GridVecMonodKokkos_PackComm
{
  int _nsubs;
  
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_float_2d _conc;
  typename AT::t_int_1d _list;
  typename AT::t_xfloat_2d _buf;
  
  GridVecMonodKokkos_PackComm(
    int nsubs,
    const typename DAT::tdual_float_2d &conc,
    const DAT::tdual_int_1d &list,
    const DAT::tdual_xfloat_1d &buf):
    _nsubs(nsubs),
    _conc(conc.view<DeviceType>()),
    _list(list.view<DeviceType>())
  {
    _buf = typename AT::t_xfloat_2d_um(
      buf.view<DeviceType>().data(),
      buf.view<DeviceType>().extent(0)/nsubs,
      nsubs);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    int cell = _list(i);
    for (int j = 0; j < _nsubs; j++) {
      _buf(i, j) = _conc(j, cell);
    }
  }
};

/* ---------------------------------------------------------------------- */

int GridVecChemostatKokkos::pack_comm_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{
  if (lmp->kokkos->forward_comm_on_host) {
    gridKK->sync(Host, CONC_MASK);
    struct GridVecMonodKokkos_PackComm<LMPHostType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPHostType>(first, last), f);
  } else {
    gridKK->sync(Device, CONC_MASK);
    struct GridVecMonodKokkos_PackComm<LMPDeviceType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType>(first, last), f);
  }
  return (last-first)*grid->nsubs;
}

/* ---------------------------------------------------------------------- */
template <class DeviceType>
struct GridVecMonodKokkos_UnpackComm
{
  int _nsubs;
  
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_float_2d _conc;
  typename AT::t_int_1d _list;
  typename AT::t_xfloat_2d _buf;
  
  GridVecMonodKokkos_UnpackComm(
    int nsubs,
    const typename DAT::tdual_float_2d &conc,
    const DAT::tdual_int_1d &list,
    const DAT::tdual_xfloat_1d &buf):
    _nsubs(nsubs),
    _conc(conc.view<DeviceType>()),
    _list(list.view<DeviceType>())
  {
    _buf = typename AT::t_xfloat_2d_um(
      buf.view<DeviceType>().data(),
      buf.view<DeviceType>().extent(0)/nsubs,
      nsubs);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    int cell = _list(i);
    for (int j = 0; j < _nsubs; j++) {
      _conc(j, cell) = _buf(i, j);
    }
  }
};

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::unpack_comm_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{
  if (lmp->kokkos->forward_comm_on_host) {
    struct GridVecMonodKokkos_UnpackComm<LMPHostType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType>(first, last), f);
    gridKK->modified(Host, CONC_MASK);
  } else {
    gridKK->sync(Device, CONC_MASK);
    struct GridVecMonodKokkos_UnpackComm<LMPDeviceType> f(grid->nsubs, gridKK->k_conc, list, buf);
    Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType>(first, last), f);
    gridKK->modified(Device, CONC_MASK);
  }
}

/* ---------------------------------------------------------------------- */

int GridVecChemostatKokkos::pack_exchange_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::unpack_exchange_kokkos(int first, int last, const DAT::tdual_int_1d &list, const DAT::tdual_xfloat_1d &buf)
{

}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::set(int narg, char **arg)
{
  sync(Host, GMASK_MASK);

  if (narg != 7) error->all(FLERR, "Invalid grid_modify set command");
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

  bulk[isub] = utils::numeric(FLERR,arg[6],true,lmp);
  if (bulk[isub] < 0) error->all(FLERR, "Illegal initial bulk concentration");

  set_grid(isub, domain, bulk[isub]);

  modified(Host, CONC_MASK | BOUNDARY_MASK | BULK_MASK | REAC_MASK);
}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::set_grid(int isub, double domain, double bulk)
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

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & GMASK_MASK) gridKK->k_mask.sync<LMPDeviceType>();
    if (mask & CONC_MASK) gridKK->k_conc.sync<LMPDeviceType>();
    if (mask & REAC_MASK) gridKK->k_reac.sync<LMPDeviceType>();
    if (mask & DENS_MASK) gridKK->k_dens.sync<LMPDeviceType>();
    if (mask & DIFF_COEFF_MASK) gridKK->k_diff_coeff.sync<LMPDeviceType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.sync<LMPDeviceType>();
    if (mask & BULK_MASK) gridKK->k_bulk.sync<LMPDeviceType>();
    if (mask & BOUNDARY_MASK) gridKK->k_boundary.sync<LMPDeviceType>();
  } else {
    if (mask & GMASK_MASK) gridKK->k_mask.sync<LMPHostType>();
    if (mask & CONC_MASK) gridKK->k_conc.sync<LMPHostType>();
    if (mask & REAC_MASK) gridKK->k_reac.sync<LMPHostType>();
    if (mask & DIFF_COEFF_MASK) gridKK->k_diff_coeff.sync<LMPHostType>();
    if (mask & DENS_MASK) gridKK->k_dens.sync<LMPHostType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.sync<LMPHostType>();
    if (mask & BULK_MASK) gridKK->k_bulk.sync<LMPHostType>();
    if (mask & BOUNDARY_MASK) gridKK->k_boundary.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (mask & GMASK_MASK) gridKK->k_mask.modify<LMPDeviceType>();
    if (mask & CONC_MASK) gridKK->k_conc.modify<LMPDeviceType>();
    if (mask & REAC_MASK) gridKK->k_reac.modify<LMPDeviceType>();
    if (mask & DIFF_COEFF_MASK) gridKK->k_diff_coeff.modify<LMPDeviceType>();
    if (mask & DENS_MASK) gridKK->k_dens.modify<LMPDeviceType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.modify<LMPDeviceType>();
    if (mask & BULK_MASK) gridKK->k_bulk.modify<LMPDeviceType>();
    if (mask & BOUNDARY_MASK) gridKK->k_boundary.modify<LMPDeviceType>();
  } else {
    if (mask & GMASK_MASK) gridKK->k_mask.modify<LMPHostType>();
    if (mask & CONC_MASK) gridKK->k_conc.modify<LMPHostType>();
    if (mask & REAC_MASK) gridKK->k_reac.modify<LMPHostType>();
    if (mask & DIFF_COEFF_MASK) gridKK->k_diff_coeff.modify<LMPHostType>();
    if (mask & DENS_MASK) gridKK->k_dens.modify<LMPHostType>();
    if (mask & GROWTH_MASK) gridKK->k_growth.modify<LMPHostType>();
    if (mask & BULK_MASK) gridKK->k_bulk.modify<LMPHostType>();
    if (mask & BOUNDARY_MASK) gridKK->k_boundary.modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void GridVecChemostatKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
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
    if ((mask & BULK_MASK) && gridKK->k_bulk.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_1d>(gridKK->k_bulk, space);
    if ((mask & BOUNDARY_MASK) && gridKK->k_boundary.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_int_2d>(gridKK->k_boundary, space);
    if ((mask & DIFF_COEFF_MASK) && gridKK->k_diff_coeff.need_sync<LMPDeviceType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_diff_coeff, space);
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
    if ((mask & BULK_MASK) && gridKK->k_bulk.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_1d>(gridKK->k_bulk, space);
    if ((mask & BOUNDARY_MASK) && gridKK->k_boundary.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_int_2d>(gridKK->k_boundary, space);
    if ((mask & DIFF_COEFF_MASK) && gridKK->k_diff_coeff.need_sync<LMPHostType>())
      perform_async_copy<DAT::tdual_float_2d>(gridKK->k_diff_coeff, space);
  }
}
