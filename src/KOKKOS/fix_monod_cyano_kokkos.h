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

FixStyle(nufeb/monod/cyano/kk,FixMonodCyanoKokkos<LMPDeviceType>)
FixStyle(nufeb/monod/cyano/kk/device,FixMonodCyanoKokkos<LMPDeviceType>)
FixStyle(nufeb/monod/cyano/kk/host,FixMonodCyanoKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_MONOD_CYANO_KOKKOS_H
#define LMP_FIX_MONOD_CYANO_KOKKOS_H

#include "fix_monod_cyano.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <int, int>
struct FixMonodCyanoCellsTag {};
struct FixMonodCyanoAtomsTag {};

template <class DeviceType>
class FixMonodCyanoKokkos: public FixMonodCyano {
 public:
  FixMonodCyanoKokkos(class LAMMPS *, int, char **);
  virtual ~FixMonodCyanoKokkos() {}
  virtual void compute();

  template <int, int> void update_cells();
  void update_atoms();

  struct Functor
  {
    int igroup;
    int ilight;   // light
    int ico2;     // co2
    int igco2;    // gas co2
    int isuc;     // sucrose
    int io2;      // dissolved co2
    
    double light_affinity;
    double co2_affinity;

    double growth;
    double yield;
    double maintain;
    double decay;
    double suc_exp;
    double gco2_flag;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_2d d_conc;
    typename AT::t_float_2d d_reac;
    typename AT::t_float_2d d_dens;
    typename AT::t_float_3d d_growth;

    Functor(FixMonodCyanoKokkos *ptr);
    
    template <int Reaction, int Growth>
    KOKKOS_INLINE_FUNCTION
    void operator()(FixMonodCyanoCellsTag<Reaction, Growth>, int) const;
  };

 protected:
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_int_1d d_mask;
  typename AT::t_float_2d d_conc;
  typename AT::t_float_2d d_reac;
  typename AT::t_float_2d d_dens;
  typename AT::t_float_3d d_growth;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
