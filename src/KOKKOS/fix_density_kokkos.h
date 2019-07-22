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

FixStyle(nufeb/density/kk,FixDensityKokkos<LMPDeviceType>)
FixStyle(nufeb/density/kk/device,FixDensityKokkos<LMPDeviceType>)
FixStyle(nufeb/density/kk/host,FixDensityKokkos<LMPHostTYpe>)

#else

#ifndef LMP_FIX_DENSITY_KOKKOS_H
#define LMP_FIX_DENSITY_KOKKOS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDensityKokkos : public Fix {
 public:
  int compute_flag;

  FixDensityKokkos(class LAMMPS *, int, char **);
  virtual ~FixDensityKokkos() {}
  int setmask();
  int modify_param(int, char **);
  virtual void post_integrate();
  virtual void compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
