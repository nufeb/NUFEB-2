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

FixStyle(nufeb/density/kk,FixDensityKokkos<LMPDeviceType>)
FixStyle(nufeb/density/kk/device,FixDensityKokkos<LMPDeviceType>)
FixStyle(nufeb/density/kk/host,FixDensityKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_DENSITY_KOKKOS_H
#define LMP_FIX_DENSITY_KOKKOS_H

#include "kokkos_type.h"
#include "fix_density.h"

namespace LAMMPS_NS {

struct FixDensityInitialTag {};
struct FixDensityComputeTag {};

template <class DeviceType>
class FixDensityKokkos : public FixDensity {
 public:
  FixDensityKokkos(class LAMMPS *, int, char **);
  ~FixDensityKokkos() {}

  void post_physics_nufeb();
  void init();
  void compute();

  struct Functor
  {
    int ngroup;
    double boxlo[3];
    double sublo[3];
    double subhi[3];
    int grid_sublo[3];
    int grid_subbox[3];
    double cell_size;
    double vol;

    typedef ArrayTypes<DeviceType> AT;

    typename AT::t_int_1d d_bitmask;
    typename AT::t_int_1d d_mask;
    typename AT::t_x_array d_x;
    typename AT::t_float_1d d_rmass;
    typename AT::t_float_1d d_biomass;
    typename AT::t_float_2d d_dens;

    Functor(FixDensityKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixDensityInitialTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixDensityComputeTag, int) const;
  };

 protected:

  double boxlo[3];
  double sublo[3];
  double subhi[3];
  int grid_sublo[3];
  int grid_subbox[3];
  double cell_size;
  double vol;
  int ngroup;

  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d d_bitmask;
  typename AT::t_int_1d d_mask;
  typename AT::t_x_array d_x;
  typename AT::t_float_1d d_rmass;
  typename AT::t_float_1d d_biomass;
  typename AT::t_float_2d d_dens;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
