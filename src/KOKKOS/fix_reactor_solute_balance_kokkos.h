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

FixStyle(nufeb/reactor/solute_balance/kk,FixReactorSoluteBalanceKokkos<LMPDeviceType>)
FixStyle(nufeb/reactor/solute_balance/kk/device,FixReactorSoluteBalanceKokkos<LMPDeviceType>)
FixStyle(nufeb/reactor/solute_balance/kk/host,FixReactorSoluteBalanceKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_REACTOR_SOLUTE_BALANCE_KOKKOS_H
#define LMP_FIX_REACTOR_SOLUTE_BALANCE_KOKKOS_H

#include "kokkos_type.h"
#include "fix_reactor_solute_balance.h"

namespace LAMMPS_NS {

struct FixSoluteBalanceComputeTag {};

template <class DeviceType>
class FixReactorSoluteBalanceKokkos : public FixReactorSoluteBalance {
 public:
  FixReactorSoluteBalanceKokkos(class LAMMPS *, int, char **);
  ~FixReactorSoluteBalanceKokkos() {}

  void compute();

  struct Functor
  {
    int iliq;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_2d d_reac;

    Functor(FixReactorSoluteBalanceKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixSoluteBalanceComputeTag, int, double &) const;
  };

 protected:
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_int_1d d_mask;
  typename AT::t_float_2d d_reac;
  HAT::t_float_1d h_bulk;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
