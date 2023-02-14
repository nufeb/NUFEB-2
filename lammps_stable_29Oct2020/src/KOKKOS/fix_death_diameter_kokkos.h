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

FixStyle(nufeb/death/diameter/kk,FixDeathDiameterKokkos<LMPDeviceType>)
FixStyle(nufeb/death/diameter/kk/device,FixDeathDiameterKokkos<LMPDeviceType>)
FixStyle(nufeb/death/diameter/kk/host,FixDeathDiameterKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_DEATH_DIAMETER_KOKKOS_H
#define LMP_FIX_DEATH_DIAMETER_KOKKOS_H

#include "kokkos_type.h"
#include "fix_death_diameter.h"

namespace LAMMPS_NS {

struct FixDeathDiameterComputeTag {};

template <class DeviceType>
class FixDeathDiameterKokkos : public FixDeathDiameter {
 public:
  FixDeathDiameterKokkos(class LAMMPS *, int, char **);
  ~FixDeathDiameterKokkos() {}

  void compute();

  struct Functor
  {
    int groupbit;
    int idead;
    int tdead;
    double diameter;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_int_1d d_type;
    typename AT::t_float_1d d_radius;

    Functor(FixDeathDiameterKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixDeathDiameterComputeTag, int) const;
  };

 private:
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d d_mask;
  typename AT::t_int_1d d_type;
  typename AT::t_float_1d d_radius;
};
}

#endif
#endif

/* ERROR/WARNING messages:
*/
