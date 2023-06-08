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

FixStyle(nufeb/division/bacillus/minicell,FixDivideBacillusMinicell)

#else

#ifndef LMP_FIX_DIVIDE_BACILLUS_MINICELL_H
#define LMP_FIX_DIVIDE_BACILLUS_MINICELL_H

#include "fix_divide.h"

namespace LAMMPS_NS {

class FixDivideBacillusMinicell : public FixDivide {
 public:
  
  FixDivideBacillusMinicell(class LAMMPS *, int, char **);
  ~FixDivideBacillusMinicell();
  void compute();
  void init();
  void *extract(const char *, int &);
  
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double *birth_length;

 private:
  int imini, type;
  double prob;
  double divlength;
  double var;
  int seed;
  double maxradius;
  int divflag, conserveflag;

  class RanPark *random;
  class AtomVecBacillus *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
