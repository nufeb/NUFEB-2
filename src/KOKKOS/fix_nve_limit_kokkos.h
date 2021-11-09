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

FixStyle(nve/limit/kk,FixNVELimitKokkos<LMPDeviceType>)
FixStyle(nve/limit/kk/device,FixNVELimitKokkos<LMPDeviceType>)
FixStyle(nve/limit/kk/host,FixNVELimitKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_NVE_LIMIT_KOKKOS_H
#define LMP_FIX_NVE_LIMIT_KOKKOS_H

#include "fix_nve_limit.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct FixNVELimitInitialTag {};
struct FixNVELimitInitialRMassTag {};
struct FixNVELimitFinalTag {};
struct FixNVELimitFinalRMassTag {};

template<class DeviceType>
class FixNVELimitKokkos : public FixNVELimit {
 public:
  FixNVELimitKokkos(class LAMMPS *, int, char **);
  ~FixNVELimitKokkos() {}
  void initial_integrate(int);
  void final_integrate();
  double compute_scalar();

  struct Functor
  {
    int groupbit;
    double dtf;
    double dtv;
    double vlimitsq;
    
    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_x_array x;
    typename AT::t_v_array v;
    typename AT::t_f_array_const f;
    typename AT::t_float_1d rmass;
    typename AT::t_float_1d mass;
    typename AT::t_int_1d type;
    typename AT::t_int_1d mask;
    
    Functor(FixNVELimitKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixNVELimitInitialTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixNVELimitInitialRMassTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixNVELimitFinalTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixNVELimitFinalRMassTag, int) const;
  };
  
 private:
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_f_array_const f;
  typename AT::t_float_1d rmass;
  typename AT::t_float_1d mass;
  typename AT::t_int_1d type;
  typename AT::t_int_1d mask;
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
