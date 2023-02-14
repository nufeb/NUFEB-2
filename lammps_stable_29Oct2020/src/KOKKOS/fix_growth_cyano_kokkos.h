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

FixStyle(nufeb/growth/cyano/kk,FixGrowthCyanoKokkos<LMPDeviceType>)
FixStyle(nufeb/growth/cyano/kk/device,FixGrowthCyanoKokkos<LMPDeviceType>)
FixStyle(nufeb/growth/cyano/kk/host,FixGrowthCyanoKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_GROWTH_CYANO_KOKKOS_H
#define LMP_FIX_GROWTH_CYANO_KOKKOS_H

#include "fix_growth_cyano.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct FixGrowthCyanoCellsReactionTag {};
struct FixGrowthCyanoCellsGrowthTag {};
struct FixGrowthCyanoAtomsTag {};

template <class DeviceType>
class FixGrowthCyanoKokkos: public FixGrowthCyano {
 public:
  FixGrowthCyanoKokkos(class LAMMPS *, int, char **);
  ~FixGrowthCyanoKokkos() {}

  void update_cells();
  void update_atoms();

  struct Functor
  {
    int igroup;
    int groupbit;
    double dt;

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


    Functor(FixGrowthCyanoKokkos *ptr);
    
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthCyanoCellsReactionTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthCyanoCellsGrowthTag, int) const;
    KOKKOS_INLINE_FUNCTION
    void operator()(FixGrowthCyanoAtomsTag, int) const;
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
