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

#include "fix_diffusion_coeff_kokkos.h"
#include "grid.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "fix_diffusion_reaction.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDiffusionCoeffKokkos<DeviceType>::FixDiffusionCoeffKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDiffusionCoeff(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDiffusionCoeffKokkos<DeviceType>::compute()
{
  d_mask = gridKK->k_mask.template view<DeviceType>();
  d_diff_coeff = gridKK->k_diff_coeff.template view<DeviceType>();
  d_dens = gridKK->k_dens.template view<DeviceType>();

  vol = grid->cell_size * grid->cell_size * grid->cell_size;

  gridKK->sync(execution_space, GMASK_MASK | DIFF_COEFF_MASK | DENS_MASK);

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDiffusionCoeffComputeTag>(0, grid->ncells), f);
  DeviceType().fence();
  copymode = 0;

  gridKK->modified(execution_space, DIFF_COEFF_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDiffusionCoeffKokkos<DeviceType>::Functor::Functor(FixDiffusionCoeffKokkos<DeviceType> *ptr):
  d_mask(ptr->d_mask), d_diff_coeff(ptr->d_diff_coeff), d_dens(ptr->d_dens),
  isub(ptr->isub), coeff_flag(ptr->coeff_flag), ratio(ptr->ratio), vol(ptr->vol),
  const_coeff(ptr->const_coeff)
{
}

/* ---------------------------------------------------------------------- */


template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDiffusionCoeffKokkos<DeviceType>::Functor::operator()(FixDiffusionCoeffComputeTag, int i) const
{
  if (!(d_mask(i) & GHOST_MASK)) {
    if (coeff_flag == 0) {
      if (d_dens(0,i) > 0) {
	d_diff_coeff(isub,i) = const_coeff * ratio;
      }
    } else if (coeff_flag == 1) {
      d_diff_coeff(isub,i) = const_coeff * (1 - (0.43 * pow(d_dens(0,i)/vol,0.92)) /
      (11.19 + 0.27 * pow(d_dens(0,i)/vol,0.99)));
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixDiffusionCoeffKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixDiffusionCoeffKokkos<LMPHostType>;
#endif
}
