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

FixStyle(nufeb/monod/nob,FixMonodNOB)

#else

#ifndef LMP_FIX_MONOD_NOB_H
#define LMP_FIX_MONOD_NOB_H

#include "fix_monod.h"

namespace LAMMPS_NS {

class FixMonodNOB: public FixMonod {
 public:
  FixMonodNOB(class LAMMPS *, int, char **);
  virtual ~FixMonodNOB() {}
  virtual void compute();

 protected:
  int io2;
  int ino2;
  int ino3;
  
  double o2_affinity;
  double no2_affinity;

  double growth;
  double yield;
  double maintain;
  double decay;
  
  template <int, int> void update_cells();
  virtual void update_atoms();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
