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

FixStyle(nufeb/growth/nob/kk,FixGrowthNOBKokkos<LMPDeviceType>)
FixStyle(nufeb/growth/nob/kk/device,FixGrowthNOBKokkos<LMPDeviceType>)
FixStyle(nufeb/growth/nob/kk/host,FixGrowthNOBKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_GROWTH_NOB_KOKKOS_H
#define LMP_FIX_GROWTH_NOB_KOKKOS_H

#include "fix_growth_nob.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct FixGrowthNOBCellsReactionTag {};
struct FixGrowthNOBCellsGrowthTag {};
struct FixGrowthNOBAtomsTag {};

template <class DeviceType>
class FixGrowthNOBKokkos: public FixGrowthNOB {
 public:
  FixGrowthNOBKokkos(class LAMMPS *, int, char **);
  ~FixGrowthNOBKokkos() {}

  void update_cells();
  void update_atoms();

  struct Functor
  {
    int igroup;
    int groupbit;
    double dt;

    int io2;
    int ino2;
    int ino3;
    
    double o2_affinity;
    double no2_affinity;

    double growth;
    double yield;
    double maintain;
    double decay;

    double boxlo[3];
    int grid_sublo[3];
    int grid_subbox[3];
    double cell_size;
    double vol;

    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_int_1d d_gmask;
    typename AT::t_float_2d d_conc;
    typename AT::t_float_2d d_reac;
    typename AT::t_float_2d d_dens;
    typename AT::t_float_3d d_growth;

    typename AT::t_x_array d_x;
    typename AT::t_float_1d d_rmass;
    typename AT::t_float_1d d_radius;
    typename AT::t_float_1d d_outer_mass;
    typename AT::t_float_1d d_outer_radius;

    Functor(FixGrowthNOBKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthNOBCellsReactionTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthNOBCellsGrowthTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthNOBAtomsTag, int) const;
  };

 protected:
  double boxlo[3];
  int grid_sublo[3];
  int grid_subbox[3];
  double cell_size;
  double vol;

  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_int_1d d_mask;
  typename AT::t_int_1d d_gmask;
  typename AT::t_float_2d d_conc;
  typename AT::t_float_2d d_reac;
  typename AT::t_float_2d d_dens;
  typename AT::t_float_3d d_growth;

  typename AT::t_x_array d_x;
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
