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

FixStyle(viscous/kk,FixViscousKokkos<LMPDeviceType>)
FixStyle(viscous/kk/device,FixViscousKokkos<LMPDeviceType>)
FixStyle(viscous/kk/host,FixViscousKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_VISCOUS_KOKKOS_H
#define LMP_FIX_VISCOUS_KOKKOS_H

#include "fix_viscous.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixViscousKokkos : public FixViscous {
 public:
  FixViscousKokkos(class LAMMPS *, int, char **);
  virtual ~FixViscousKokkos() {}
  void post_force(int);

  struct Functor
  {
    int groupbit;
    
    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_type;
    typename AT::t_int_1d d_mask;
    typename AT::t_v_array d_v;
    typename AT::t_f_array d_f;
    typename AT::t_float_1d d_gamma;

    Functor(FixViscousKokkos *ptr);
    
    KOKKOS_INLINE_FUNCTION
    void operator()(int i) const;
  };
  
 protected:
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d d_type;
  typename AT::t_int_1d d_mask;
  typename AT::t_v_array d_v;
  typename AT::t_f_array d_f;
  typename AT::t_float_1d d_gamma;
  typename HAT::t_float_1d h_gamma;
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
