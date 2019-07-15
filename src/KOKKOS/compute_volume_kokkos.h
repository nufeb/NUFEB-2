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

#ifdef COMPUTE_CLASS

ComputeStyle(nufeb/volume/kk,ComputeVolumeKokkos<LMPDeviceType>)
ComputeStyle(nufeb/volume/kk/device,ComputeVolumeKokkos<LMPDeviceType>)
ComputeStyle(nufeb/volume/kk/host,ComputeVolumeKokkos<LMPHostType>)

#else

#ifndef LMP_COMPUTE_VOLUME_KOKKOS_H
#define LMP_COMPUTE_VOLUME_KOKKOS_H

#include "compute_volume.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class ComputeVolumeKokkos : public ComputeVolume {
 public:
  ComputeVolumeKokkos(class LAMMPS *, int, char **);
  virtual ~ComputeVolumeKokkos() {}
  virtual double compute_scalar();

  struct Functor
  {
    int groupbit;
    
    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_1d d_radius;
    
    Functor(ComputeVolumeKokkos *ptr);
    
    KOKKOS_INLINE_FUNCTION
    void operator()(int, double &) const;
  };
  
 protected:
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d d_mask;
  typename AT::t_float_1d d_radius;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
