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

FixStyle(nufeb/density,FixDensity)

#else

#ifndef LMP_FIX_DENSITY_H
#define LMP_FIX_DENSITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDensity : public Fix {
 public:
  FixDensity(class LAMMPS *, int, char **);
  virtual ~FixDensity() {}
  int setmask();
  int modify_param(int, char **);
  virtual void post_physics_nufeb();
  virtual void compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
