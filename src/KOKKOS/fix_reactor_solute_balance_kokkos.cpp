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

#include "fix_reactor_solute_balance_kokkos.h"
#include "grid.h"
#include "update.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixReactorSoluteBalanceKokkos<DeviceType>::FixReactorSoluteBalanceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixReactorSoluteBalance(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixReactorSoluteBalanceKokkos<DeviceType>::compute()
{
  double vol = grid->cell_size * grid->cell_size * grid->cell_size;

  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_reac = gridKK->k_reac.template view<DeviceType>();
  h_bulk = gridKK->k_bulk.view<LMPHostType>();

  gridKK->sync(execution_space, GMASK_MASK | REAC_MASK);

  copymode = 1;
  Functor f(this);
  double result = 0.0;
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<
    DeviceType,
    FixSoluteBalanceComputeTag>(0, grid->ncells), f, Kokkos::Sum<double>(result));
  DeviceType().fence();
  copymode = 0;

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);

  h_bulk(iliq) += ((q / rvol) * (inlet - h_bulk(iliq)) +
      ((reactor_af * result * vol) / (rvol * domain_af))) * update->dt;
  h_bulk(iliq) = MAX(0, h_bulk(iliq));
  gridKK->modified(Host, BULK_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixReactorSoluteBalanceKokkos<DeviceType>::Functor::Functor(FixReactorSoluteBalanceKokkos<DeviceType> *ptr):
  d_mask(ptr->d_mask), d_reac(ptr->d_reac), iliq(ptr->iliq)
{}

/* ---------------------------------------------------------------------- */


template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixReactorSoluteBalanceKokkos<DeviceType>::Functor::operator()(FixSoluteBalanceComputeTag, int i, double &sum) const
{
  if (!(d_mask(i) & GHOST_MASK)) {
    sum += d_reac(iliq, i);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixReactorSoluteBalanceKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixReactorSoluteBalanceKokkos<LMPHostType>;
#endif
}
