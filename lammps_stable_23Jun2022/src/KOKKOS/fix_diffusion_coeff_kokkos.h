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

FixStyle(nufeb/diffusion_coeff/kk,FixDiffusionCoeffKokkos<LMPDeviceType>)
FixStyle(nufeb/diffusion_coeff/kk/device,FixDiffusionCoeffKokkos<LMPDeviceType>)
FixStyle(nufeb/diffusion_coeff/kk/host,FixDiffusionCoeffKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_DIFFUSION_COEFF_KOKKOS_H
#define LMP_FIX_DIFFUSION_COEFF_KOKKOS_H

#include "kokkos_type.h"
#include "fix_diffusion_coeff.h"

namespace LAMMPS_NS {

struct FixDiffusionCoeffComputeTag {};

template <class DeviceType>
class FixDiffusionCoeffKokkos : public FixDiffusionCoeff {
 public:
  FixDiffusionCoeffKokkos(class LAMMPS *, int, char **);
  ~FixDiffusionCoeffKokkos() {}

  void compute();

  struct Functor
  {
    int isub;
    int coeff_flag;
    double ratio;
    double vol;
    double const_coeff;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_2d d_diff_coeff;
    typename AT::t_float_2d d_dens;

    Functor(FixDiffusionCoeffKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixDiffusionCoeffComputeTag, int) const;
  };

 protected:
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d d_mask;
  typename AT::t_float_2d d_diff_coeff;
  typename AT::t_float_2d d_dens;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
