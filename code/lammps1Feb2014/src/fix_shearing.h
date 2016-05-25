/* ----------------------------------------------------------------------
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

FixStyle(shearing,FixShearing)

#else

#ifndef LMP_FIX_SHEARING_H
#define LMP_FIX_SHEARING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixShearing : public Fix {
 public:
	FixShearing(class LAMMPS *, int, char **);
 ~FixShearing();
  int setmask();
  void init();
  virtual void post_force(int);

 private:
  char **var;
  int *ivar;
  bigint tmin, tmax;
  double visco;
  double rate;
  double dflag;
};

}

#endif
#endif


