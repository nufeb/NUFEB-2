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

FixStyle(nufeb/growth/anammox,FixGrowthAnammox)

#else

#ifndef LMP_FIX_GROWTH_ANAMMOX_H
#define LMP_FIX_GROWTH_ANAMMOX_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthAnammox: public FixGrowth {
 public:
  FixGrowthAnammox(class LAMMPS *, int, char **);
  virtual ~FixGrowthAnammox() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  int inh4;
  int io2;
  int ino2;
  int ino3;
  int in2;
  
  double nh4_affinity;
  double o2_affinity;
  double no2_affinity;

  double growth;
  double yield;
  double maintain;
  double decay;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
