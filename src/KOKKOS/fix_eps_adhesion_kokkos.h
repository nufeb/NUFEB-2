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

FixStyle(nufeb/adhesion/eps/kk,FixEPSAdhesionKokkos<LMPDeviceType>)
FixStyle(nufeb/adhesion/eps/kk/device,FixEPSAdhesionKokkos<LMPDeviceType>)
FixStyle(nufeb/adhesion/eps/kk/host,FixEPSAdhesionKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_EPS_ADHESION_KOKKOS_H
#define LMP_FIX_EPS_ADHESION_KOKKOS_H

#include "fix_adhesion_eps.h"
#include "kokkos_type.h"
  
namespace LAMMPS_NS {

template <int, int, int>
struct FixEPSAdhesionTag {};
  
template <class DeviceType>
class FixEPSAdhesionKokkos : public FixEPSAdhesion {
 public:
  FixEPSAdhesionKokkos(class LAMMPS *, int, char **);
  ~FixEPSAdhesionKokkos() {}
  virtual void post_force(int);

  struct Functor
  {
    int groupbit;
    int epsmask;
    int nlocal;
    double ke;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_neighbors_2d d_neighbors;
    typename AT::t_int_1d_randomread d_ilist;
    typename AT::t_int_1d_randomread d_numneigh;
    typename AT::t_x_array d_x;
    typename AT::t_f_array d_f;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_1d d_rmass;
    typename AT::t_float_1d d_radius;
    typename AT::t_float_1d d_outer_mass;
    typename AT::t_float_1d d_outer_radius;

    Functor(FixEPSAdhesionKokkos *ptr);

    template <int NEIGHFLAG, int NEWTON_PAIR, int DISP>
    KOKKOS_INLINE_FUNCTION
    void operator()(FixEPSAdhesionTag<NEIGHFLAG, NEWTON_PAIR, DISP>, int) const;
  };

 private:
  int epsmask;
  int nlocal;
  
  typedef ArrayTypes<DeviceType> AT;

  // for neighbor list lookup
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_x_array d_x;
  typename AT::t_f_array d_f;
  typename AT::t_int_1d d_mask;
  typename AT::t_float_1d d_rmass;
  typename AT::t_float_1d d_radius;
  typename AT::t_float_1d d_outer_mass;
  typename AT::t_float_1d d_outer_radius;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
