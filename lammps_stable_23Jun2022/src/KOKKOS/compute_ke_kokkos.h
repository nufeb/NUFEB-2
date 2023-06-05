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

#ifdef COMPUTE_CLASS

ComputeStyle(ke/kk,ComputeKEKokkos<LMPDeviceType>)
ComputeStyle(ke/kk/device,ComputeKEKokkos<LMPDeviceType>)
ComputeStyle(ke/kk/host,ComputeKEKokkos<LMPHostType>)

#else

#ifndef LMP_COMPUTE_KE_KOKKOS_H
#define LMP_COMPUTE_KE_KOKKOS_H

#include "compute_ke.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct ComputeKEMassTag {};
struct ComputeKERMassTag {};
  
template <class DeviceType>
class ComputeKEKokkos : public ComputeKE {
 public:
  ComputeKEKokkos(class LAMMPS *, int, char **);
  double compute_scalar();

  struct Functor
  {
    int groupbit;
    
    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_type;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_1d d_rmass;
    typename AT::t_float_1d d_mass;
    typename AT::t_v_array d_v;

    Functor(ComputeKEKokkos *ptr);
    
    KOKKOS_INLINE_FUNCTION
    void operator()(ComputeKEMassTag, int, double &) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(ComputeKERMassTag, int, double &) const;
  };
    
 protected:
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d d_type;
  typename AT::t_int_1d d_mask;
  typename AT::t_float_1d d_rmass;
  typename AT::t_float_1d d_mass;
  typename AT::t_v_array d_v;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
