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

#ifdef FIX_CLASS

FixStyle(nufeb/division/coccus/kk,FixDivideCoccusKokkos<LMPDeviceType>)
FixStyle(nufeb/division/coccus/kk/device,FixDivideCoccusKokkos<LMPDeviceType>)
FixStyle(nufeb/division/coccus/kk/host,FixDivideCoccusKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_DIVIDE_COCCUS_KOKKOS_H
#define LMP_FIX_DIVIDE_COCCUS_KOKKOS_H

#include "fix_divide_coccus.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"
#include "comm_kokkos.h"
  
namespace LAMMPS_NS {

struct FixDivideCoccusComputeTag {};
  
template <class DeviceType>
class FixDivideCoccusKokkos : public FixDivideCoccus {
 public:
  FixDivideCoccusKokkos(class LAMMPS *, int, char **);
  ~FixDivideCoccusKokkos();
  void init();

  void compute();
  void grow_arrays(int);

  struct Functor
  {
    int groupbit;
    double eps_density;

    double boxlo[3];
    double boxhi[3];

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_divide_list;

    typename AT::t_x_array d_x;
    typename AT::t_f_array d_f;
    typename AT::t_v_array d_v;
    typename AT::t_tagint_1d d_tag;
    typename AT::t_int_1d d_mask;
    typename AT::t_v_array d_omega;
    typename AT::t_f_array d_torque;

    typename AT::t_float_1d d_rmass;
    typename AT::t_float_1d d_biomass;
    typename AT::t_float_1d d_radius;
    typename AT::t_float_1d d_outer_mass;
    typename AT::t_float_1d d_outer_radius;

    Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
    typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;

    Functor(FixDivideCoccusKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixDivideCoccusComputeTag, int) const;
  };

 private:
  int nlocal;
  double boxlo[3];
  double boxhi[3];
  int *divide_list;
  
  typedef ArrayTypes<DeviceType> AT;
  typename AT::tdual_int_1d k_divide_list;
  typename AT::t_int_1d d_divide_list;

  typename AT::t_x_array d_x;
  typename AT::t_f_array d_f;
  typename AT::t_v_array d_v;
  typename AT::t_tagint_1d d_tag;
  typename AT::t_int_1d d_mask;
  typename AT::t_v_array d_omega;
  typename AT::t_f_array d_torque;

  typename AT::t_float_1d d_rmass;
  typename AT::t_float_1d d_biomass;
  typename AT::t_float_1d d_radius;
  typename AT::t_float_1d d_outer_mass;
  typename AT::t_float_1d d_outer_radius;

  HAT::t_int_1d h_divide_list;
  HAT::t_int_1d h_mask;
  HAT::t_float_1d h_radius;
  HAT::t_int_1d h_type;

  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
