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

#include <mpi.h>
#include "compute_ke_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeKEKokkos<DeviceType>::ComputeKEKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeKE(lmp, narg, arg)
{
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = TYPE_MASK | MASK_MASK | RMASS_MASK | V_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double ComputeKEKokkos<DeviceType>::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  d_type = atomKK->k_type.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_mass = atomKK->k_mass.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();

  double ke = 0.0;

  copymode = 1;
  Functor f(this);
  if (atom->rmass)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, ComputeKERMassTag>(0, atom->nlocal), f, ke);
  else
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, ComputeKEMassTag>(0, atom->nlocal), f, ke);
  Kokkos::fence();
  copymode = 0;

  MPI_Allreduce(&ke,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar *= pfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeKEKokkos<DeviceType>::Functor::Functor(ComputeKEKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit),
  d_type(ptr->d_type), d_mask(ptr->d_mask), d_rmass(ptr->d_rmass),
  d_mass(ptr->d_mass), d_v(ptr->d_v) {}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeKEKokkos<DeviceType>::Functor::operator()(ComputeKEMassTag, int i, double &ke) const
{
  if (d_mask[i] & groupbit)
    ke += d_mass[d_type[i]] *
      (d_v(i, 0)*d_v(i, 0) + d_v(i, 1)*d_v(i, 1) + d_v(i, 2)*d_v(i, 2));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeKEKokkos<DeviceType>::Functor::operator()(ComputeKERMassTag, int i, double &ke) const
{
  if (d_mask[i] & groupbit)
    ke += d_rmass[i] * (d_v(i, 0)*d_v(i, 0) + d_v(i, 1)*d_v(i, 1) + d_v(i, 2)*d_v(i, 2));
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeKEKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class ComputeKEKokkos<LMPHostType>;
#endif
}
