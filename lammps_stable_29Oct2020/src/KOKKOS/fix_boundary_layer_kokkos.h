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

#ifdef FIX_CLASS

FixStyle(nufeb/boundary_layer/kk,FixBoundaryLayerKokkos<LMPDeviceType>)
FixStyle(nufeb/boundary_layer/kk/device,FixBoundaryLayerKokkos<LMPDeviceType>)
FixStyle(nufeb/boundary_layer/kk/host,FixBoundaryLayerKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_BOUNDARY_LAYER_KOKKOS_H
#define LMP_FIX_BOUNDARY_LAYER_KOKKOS_H

#include "kokkos_type.h"
#include "fix_boundary_layer.h"

namespace LAMMPS_NS {

struct FixLayerUpdateMaskTag {};
struct FixLayerComputeMinxTag {};
struct FixLayerComputeMinyTag {};
struct FixLayerComputeMinzTag {};
struct FixLayerComputeMaxxTag {};
struct FixLayerComputeMaxyTag {};
struct FixLayerComputeMaxzTag {};

template <class DeviceType>
class FixBoundaryLayerKokkos : public FixBoundaryLayer {
 public:
  FixBoundaryLayerKokkos(class LAMMPS *, int, char **);
  ~FixBoundaryLayerKokkos() {}

  void compute();

  struct Functor
  {
    int isub;
    int groupbit;

    int grid_box[3];
    int grid_subbox[3];
    double grid_sublo[3];
    double grid_subhi[3];

    int layerhi[3];
    int layerlo[3];
    int sublayerhi[3];
    int sublayerlo[3];

    typedef ArrayTypes<DeviceType> AT;

    typename AT::t_x_array d_x;
    typename AT::t_int_1d d_mask;
    typename AT::t_int_1d d_gmask;
    typename AT::t_int_2d d_boundary;

    Functor(FixBoundaryLayerKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerComputeMinxTag, int, double &) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerComputeMinyTag, int, double &) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerComputeMinzTag, int, double &) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerComputeMaxxTag, int, double &) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerComputeMaxyTag, int, double &) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerComputeMaxzTag, int, double &) const;

    KOKKOS_INLINE_FUNCTION
    void operator()(FixLayerUpdateMaskTag, int) const;
  };

 protected:
  int grid_box[3];
  int grid_subbox[3];
  double grid_sublo[3];
  double grid_subhi[3];

  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_x_array d_x;
  typename AT::t_int_1d d_mask;
  typename AT::t_int_1d d_gmask;
  typename AT::t_int_2d d_boundary;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
