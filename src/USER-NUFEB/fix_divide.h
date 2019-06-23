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

FixStyle(nufeb/divide,FixDivide)

#else

#ifndef LMP_FIX_DIVIDE_H
#define LMP_FIX_DIVIDE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDivide : public Fix {
 public:
  int compute_flag;
  
  FixDivide(class LAMMPS *, int, char **);
  ~FixDivide();
  int setmask();
  int modify_param(int, char **);
  void post_integrate();
  void post_neighbor();
  void compute();
  
 private:
  double diameter;
  double eps_density;
  int seed;  

  class RanPark *random;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
