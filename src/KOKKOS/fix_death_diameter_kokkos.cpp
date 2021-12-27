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

#include "fix_death_diameter_kokkos.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "update.h"
#include "atom_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixDeathDiameterKokkos<DeviceType>::FixDeathDiameterKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDeathDiameter(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = MASK_MASK | RADIUS_MASK | TYPE_MASK;
  datamask_modify = MASK_MASK | TYPE_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixDeathDiameterKokkos<DeviceType>::compute()
{
  atomKK->sync(execution_space,datamask_read);

  d_mask = atomKK->k_mask.view<DeviceType>();
  d_type = atomKK->k_type.view<DeviceType>();
  d_radius = atomKK->k_radius.template view<DeviceType>();

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDeathDiameterComputeTag>(0, atomKK->nlocal), f);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixDeathDiameterKokkos<DeviceType>::Functor::Functor(FixDeathDiameterKokkos<DeviceType> *ptr):
  d_mask(ptr->d_mask), groupbit(ptr->groupbit), idead(ptr->idead),
  diameter(ptr->diameter), d_radius(ptr->d_radius), tdead(ptr->tdead),
  d_type(ptr->d_type)
{}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDeathDiameterKokkos<DeviceType>::Functor::operator()(FixDeathDiameterComputeTag, int i) const
{
  if (d_mask(i) & groupbit) {
    if (d_radius(i) < 0.5 * diameter) {
      d_mask(i) = idead;
      if (tdead > 0) d_type(i) = tdead;
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixDeathDiameterKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixDeathDiameterKokkos<LMPHostType>;
#endif
}
