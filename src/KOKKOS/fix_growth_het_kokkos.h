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

FixStyle(nufeb/monod/het/kk,FixMonodHETKokkos<LMPDeviceType>)
FixStyle(nufeb/monod/het/kk/device,FixMonodHETKokkos<LMPDeviceType>)
FixStyle(nufeb/monod/het/kk/host,FixMonodHETKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_MONOD_HET_KOKKOS_H
#define LMP_FIX_MONOD_HET_KOKKOS_H

#include "../USER-NUFEB/fix_growth_het.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <int, int>
struct FixMonodHETCellsTag {};
struct FixMonodHETAtomsTag {};

template <class DeviceType>
class FixMonodHETKokkos: public FixMonodHET {
 public:
  FixMonodHETKokkos(class LAMMPS *, int, char **);
  virtual ~FixMonodHETKokkos() {}
  virtual void compute();

  template <int, int> void update_cells();
  void update_atoms();

  struct Functor
  {
    int igroup;
    int isub;
    int io2;
    int ino2;
    int ino3;
    
    double sub_affinity;
    double o2_affinity;
    double no2_affinity;
    double no3_affinity;

    double growth;
    double yield;
    double maintain;
    double decay;
    double eps_yield;
    double anoxic;
    double eps_dens;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_2d d_conc;
    typename AT::t_float_2d d_reac;
    typename AT::t_float_2d d_dens;
    typename AT::t_float_3d d_growth;

    Functor(FixMonodHETKokkos *ptr);
    
    template <int Reaction, int Growth>
    KOKKOS_INLINE_FUNCTION
    void operator()(FixMonodHETCellsTag<Reaction, Growth>, int) const;
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
