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

FixStyle(nufeb/monod/simple,FixMonodSimple)

#else

#ifndef LMP_FIX_MONOD_SIMPLE_H
#define LMP_FIX_MONOD_SIMPLE_H

#include "fix_monod.h"

namespace LAMMPS_NS {

class FixMonodSimple: public FixMonod {
 public:
  FixMonodSimple(class LAMMPS *, int, char **);
  virtual ~FixMonodSimple() {}
  virtual void compute();

 protected:
  int isub;
  double sub_affinity;

  double growth;
  double yield;
  double maintain;
  double decay;
  
  class AtomVecBacillus *avec;

  template <int, int> void update_cells();
  virtual void update_atoms();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
