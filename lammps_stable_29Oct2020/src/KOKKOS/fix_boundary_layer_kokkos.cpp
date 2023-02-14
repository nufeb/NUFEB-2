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

#include "fix_boundary_layer_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "domain.h"
#include "grid.h"
#include "grid_kokkos.h"
#include "grid_masks.h"
#include "group.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixBoundaryLayerKokkos<DeviceType>::FixBoundaryLayerKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixBoundaryLayer(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | MASK_MASK;
  datamask_modify = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixBoundaryLayerKokkos<DeviceType>::compute()
{
  d_mask = atomKK->k_mask.template view<DeviceType>();
  d_x = atomKK->k_x.template view<DeviceType>();
  d_boundary = gridKK->k_boundary.template view<DeviceType>();
  d_gmask = gridKK->k_mask.template view<DeviceType>();

  atomKK->sync(execution_space, datamask_read);
  gridKK->sync(execution_space, GMASK_MASK | BOUNDARY_MASK);

  copymode = 1;
  Functor f(this);
  double extremum[6];
  extremum[0] = domain->prd[0];
  extremum[1] = domain->prd[1];
  extremum[2] = domain->prd[2];
  extremum[3] = extremum[4] = extremum[5] = 0.0;

  // compute maximum and minimum atom position w.r.t 6 planes in the system
  // not use LAMMPS_LAMBDA due to incompatiability issue with some gpu devices
  // xlo
  if (boundary[0] == DIRICHLET) {
    Kokkos::parallel_reduce(
	Kokkos::RangePolicy<
	DeviceType,
	FixLayerComputeMinxTag>(0, atomKK->nlocal), f, Kokkos::Min<double>(extremum[0]));
  }

  // ylo
  if (boundary[2] == DIRICHLET) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<
      DeviceType,
      FixLayerComputeMinyTag>(0, atomKK->nlocal), f, Kokkos::Min<double>(extremum[1]));
  }

  // zlo
  if (boundary[4] == DIRICHLET) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<
      DeviceType,
      FixLayerComputeMinzTag>(0, atomKK->nlocal), f, Kokkos::Min<double>(extremum[2]));
  }

  // xhi
  if (boundary[1] == DIRICHLET) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<
      DeviceType,
      FixLayerComputeMaxxTag>(0, atomKK->nlocal), f, Kokkos::Max<double>(extremum[3]));
  }

  // yhi
  if (boundary[3] == DIRICHLET) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<
      DeviceType,
      FixLayerComputeMaxyTag>(0, atomKK->nlocal), f, Kokkos::Max<double>(extremum[4]));
  }

  // zhi
  if (boundary[5] == DIRICHLET) {
    Kokkos::parallel_reduce(
      Kokkos::RangePolicy<
      DeviceType,
      FixLayerComputeMaxzTag>(0, atomKK->nlocal), f, Kokkos::Max<double>(extremum[5]));
  }

  DeviceType().fence();
  MPI_Allreduce(MPI_IN_PLACE, &extremum[0], 3, MPI_DOUBLE, MPI_MIN, world);
  MPI_Allreduce(MPI_IN_PLACE, &extremum[3], 3, MPI_DOUBLE, MPI_MAX, world);

  for (int i = 0; i < 3; i++) {
    layerlo[i] = static_cast<int>((extremum[i] - height) / grid->cell_size);
    layerhi[i] = static_cast<int>((height + extremum[3+i]) / grid->cell_size) + 1;
    sublayerlo[i] = layerlo[i] - (grid->sublo[i] + 1);
    sublayerhi[i] = layerhi[i] - (grid->sublo[i] + 1);

    grid_box[i] = grid->box[i];
    grid_subbox[i] = grid->subbox[i];
    grid_sublo[i] = grid->sublo[i];
    grid_subhi[i] = grid->subhi[i];
  }

  Functor f1(this);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<
      DeviceType,
      FixLayerUpdateMaskTag>(0, grid->subbox[2]),f1);
  copymode = 0;

  gridKK->modified(execution_space, GMASK_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixBoundaryLayerKokkos<DeviceType>::Functor::Functor(FixBoundaryLayerKokkos<DeviceType> *ptr):
  d_mask(ptr->d_mask), d_x(ptr->d_x), d_gmask(ptr->d_gmask),
  d_boundary(ptr->d_boundary), groupbit(ptr->groupbit),
  isub(ptr->isub)
{
  for (int i = 0; i < 3; i++) {
    layerhi[i] = ptr->layerhi[i];
    layerlo[i] = ptr->layerlo[i];
    sublayerhi[i] = ptr->sublayerhi[i];
    sublayerlo[i] = ptr->sublayerlo[i];

    grid_box[i] = ptr->grid_box[i];
    grid_subbox[i] = ptr->grid_subbox[i];
    grid_sublo[i] = ptr->grid_sublo[i];
    grid_subhi[i] = ptr->grid_subhi[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerComputeMinxTag, int i, double& minx) const
{
  if (d_mask(i) & groupbit)
    minx = MIN(minx, d_x(i,0));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerComputeMinyTag, int i, double& miny) const
{
  if (d_mask(i) & groupbit)
    miny = MIN(miny, d_x(i,1));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerComputeMinzTag, int i, double& minz) const
{
  if (d_mask(i) & groupbit)
    minz = MIN(minz, d_x(i,2));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerComputeMaxxTag, int i, double& maxx) const
{
  if (d_mask(i) & groupbit)
    maxx = MAX(maxx, d_x(i,0));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerComputeMaxyTag, int i, double& maxy) const
{
  if (d_mask(i) & groupbit)
    maxy = MAX(maxy, d_x(i,1));
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerComputeMaxzTag, int i, double& maxz) const
{
  if (d_mask[i] & groupbit)
    maxz = MAX(maxz, d_x(i,2));
}
/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixBoundaryLayerKokkos<DeviceType>::Functor::operator()(FixLayerUpdateMaskTag, int z) const
{
  for (int y = 0; y < grid_subbox[1]; y++) {
    for (int x = 0; x < grid_subbox[0]; x++) {
      int i = x + y * grid_subbox[0] + z * grid_subbox[0] * grid_subbox[1];
      int m = 0;
      if (d_gmask(i) & GHOST_MASK)
	m |= GHOST_MASK;
      if (d_gmask(i) & CORNER_MASK)
	m |= CORNER_MASK;
      else {
	if ((grid_sublo[0] < 0 && x == 0) ||
	    (d_boundary(isub, 0) == 0 && x <= sublayerlo[0]))
	  m |= X_NB_MASK;
	if ((grid_subhi[0] > grid_box[0] && x == grid_subbox[0] - 1) ||
	    (d_boundary(isub, 1) == 0 && x >= sublayerhi[0]))
	  m |= X_PB_MASK;
	if ((grid_sublo[1] < 0 && y == 0) ||
	    (d_boundary(isub, 2) == 0 && y <= sublayerlo[1]))
	  m |= Y_NB_MASK;
	if ((grid_subhi[1] > grid_box[1] && y == grid_subbox[1] - 1) ||
	    (d_boundary(isub, 3) == 0 && y >= sublayerhi[1]))
	  m |= Y_PB_MASK;
	if ((grid_sublo[2] < 0 && z == 0) ||
	    (d_boundary(isub, 4) == 0 && z <= sublayerlo[2]))
	  m |= Z_NB_MASK;
	if ((grid_subhi[2] > grid_box[2] && z == grid_subbox[2] - 1) ||
	    (d_boundary(isub, 5) == 0 && z >= sublayerhi[2]))
	  m |= Z_PB_MASK;
      }
      d_gmask(i) = m;
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixBoundaryLayerKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixBoundaryLayerKokkos<LMPHostType>;
#endif
}
