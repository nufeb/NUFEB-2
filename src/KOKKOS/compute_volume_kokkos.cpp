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

#include "compute_volume_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

#define FOURTHIRDSPI 4.18879020478

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeVolumeKokkos<DeviceType>::ComputeVolumeKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeVolume(lmp, narg, arg)
{
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = RADIUS_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double ComputeVolumeKokkos<DeviceType>::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  d_mask = atomKK->k_mask.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();

  copymode = 1;
  Functor f(this);
  scalar = 0.0;
  Kokkos::parallel_reduce(atom->nlocal, f, scalar);
  Kokkos::fence();
  copymode = 0;

  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, MPI_DOUBLE, MPI_SUM, world); 
  return scalar;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeVolumeKokkos<DeviceType>::Functor::Functor(ComputeVolumeKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit),
  d_mask(ptr->d_mask), d_radius(ptr->d_radius) {}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeVolumeKokkos<DeviceType>::Functor::operator()(int i, double &sum) const
{
  if (d_mask[i] & groupbit) {
    double r = d_radius[i];
    sum += FOURTHIRDSPI * r * r * r;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeVolumeKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class ComputeVolumeKokkos<LMPHostType>;
#endif
}
