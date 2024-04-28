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

FixStyle(mutation/mutate,FixMutate)

#else

#ifndef LMP_FIX_MUTATE_H
#define LMP_FIX_MUTATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMutate : public Fix {
 public:
  FixMutate(class LAMMPS *, int, char **);
  ~FixMutate() {}
  int modify_param(int, char **);

  void init() {}
  int setmask();
  void biology_nufeb();
  void compute();

 private:
  class RanPark *random;

  int imutant;  // group index of target species
  double prob;  // mutation probability
  int seed;     // random seed
};

}

#endif
#endif
