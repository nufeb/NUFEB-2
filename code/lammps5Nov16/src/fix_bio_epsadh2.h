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

FixStyle(epsadh2,FixEPSAdh2)

#else

#ifndef LMP_FIX_EPSADH2_H
#define LMP_FIX_EPSADH2_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEPSAdh2 : public Fix {

 public:
  FixEPSAdh2(class LAMMPS *, int, char **);
 ~FixEPSAdh2();
  int setmask();
  void init();
  virtual void post_force(int);

  double smax; //maximum seperatioin for force cutoff

 private:
  char *var;
  int ivar;
  int nmax, nvalues;
  int npairs;
  bigint laststep;
  bigint laststep_local;
  int flag;

  class AtomVecBio *avec;
};

}

#endif
#endif


