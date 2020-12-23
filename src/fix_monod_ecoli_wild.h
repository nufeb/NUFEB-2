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

FixStyle(nufeb/monod/ecoli/wild,FixMonodEcoliWild)

#else

#ifndef LMP_FIX_MONOD_ECOLI_WILD_H
#define LMP_FIX_MONOD_ECOLI_WILD_H

#include "fix_monod.h"

namespace LAMMPS_NS {

class FixMonodEcoliWild: public FixMonod {
 public:
  FixMonodEcoliWild(class LAMMPS *, int, char **);
  virtual ~FixMonodEcoliWild() {}
  virtual void compute();

 protected:
  int isuc;   // sucrose
  int io2;
  int ico2;

  double suc_affinity;
  double o2_affinity;

  double growth;
  double yield;
  double maintain;
  double decay;
  
  template <int, int> void update_cells();
  void update_atoms();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
