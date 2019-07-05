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

FixStyle(nufeb/monod/het,FixMonodHET)

#else

#ifndef LMP_FIX_MONOD_HET_H
#define LMP_FIX_MONOD_HET_H

#include "fix_monod.h"

namespace LAMMPS_NS {

class FixMonodHET: public FixMonod {
 public:
  FixMonodHET(class LAMMPS *, int, char **);
  virtual ~FixMonodHET() {}
  virtual void compute();

  template <int, int> void update_cells();
  void update_atoms();

 protected:
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
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
