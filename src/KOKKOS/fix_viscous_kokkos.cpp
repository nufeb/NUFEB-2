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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_viscous_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixViscousKokkos<DeviceType>::FixViscousKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixViscous(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = TYPE_MASK | MASK_MASK | V_MASK | F_MASK;
  datamask_modify = F_MASK;

  d_gamma = typename AT::t_float_1d("viscous:d_gamma", atom->ntypes+1);
  h_gamma = typename HAT::t_float_1d("viscous:h_gamma", atom->ntypes+1);
  for (int i = 0; i <= atom->ntypes; i++)
    h_gamma[i] = gamma[i];
  Kokkos::deep_copy(d_gamma, h_gamma);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixViscousKokkos<DeviceType>::post_force(int /*vflag*/)
{
  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();

  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(atom->nlocal, f);
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixViscousKokkos<DeviceType>::Functor::Functor(FixViscousKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit),
  d_type(ptr->d_type), d_mask(ptr->d_mask),
  d_v(ptr->d_v), d_f(ptr->d_f), d_gamma(ptr->d_gamma) {}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixViscousKokkos<DeviceType>::Functor::operator()(int i) const
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type
  if (d_mask[i] & groupbit) {
    double drag = d_gamma[d_type[i]];
    d_f(i, 0) -= drag*d_v(i, 0);
    d_f(i, 1) -= drag*d_v(i, 1);
    d_f(i, 2) -= drag*d_v(i, 2);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixViscousKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixViscousKokkos<LMPHostType>;
#endif
}
