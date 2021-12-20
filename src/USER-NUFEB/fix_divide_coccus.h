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

FixStyle(nufeb/division/coccus,FixDivideCoccus)

#else

#ifndef LMP_FIX_DIVIDE_COCCUS_H
#define LMP_FIX_DIVIDE_COCCUS_H

#include "fix_divide.h"

namespace LAMMPS_NS {

class FixDivideCoccus : public FixDivide {
 public:
  
  FixDivideCoccus(class LAMMPS *, int, char **);
  virtual ~FixDivideCoccus();
  virtual void compute();
  
 protected:
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
