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

FixStyle(nufeb/growth/eps/kk,FixGrowthEPSKokkos<LMPDeviceType>)
FixStyle(nufeb/growth/eps/kk/device,FixGrowthEPSKokkos<LMPDeviceType>)
FixStyle(nufeb/growth/eps/kk/host,FixGrowthEPSKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_GROWTH_EPS_KOKKOS_H
#define LMP_FIX_GROWTH_EPS_KOKKOS_H

#include "fix_growth_eps.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct FixGrowthEPSCellsReactionTag {};
struct FixGrowthEPSCellsGrowthTag {};
struct FixGrowthEPSAtomsTag {};

template <class DeviceType>
class FixGrowthEPSKokkos: public FixGrowthEPS {
 public:
  FixGrowthEPSKokkos(class LAMMPS *, int, char **);
  ~FixGrowthEPSKokkos() {}

  void update_cells();
  void update_atoms();

  struct Functor
  {
    int igroup;
    int groupbit;
    double dt;

    int isub;
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

    Functor(FixGrowthEPSKokkos *ptr);

    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthEPSCellsReactionTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthEPSCellsGrowthTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthEPSAtomsTag, int) const;
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
