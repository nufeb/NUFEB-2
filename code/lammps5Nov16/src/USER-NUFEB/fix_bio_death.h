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

FixStyle(death,FixDeath)

#else

#ifndef LMP_FIX_DEATH_H
#define LMP_FIX_DEATH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeath : public Fix {
 public:
	FixDeath(class LAMMPS *, int, char **);
  ~FixDeath();
  int setmask();
  void init();
  void pre_exchange();

 private:
  class AtomVecBio *avec;

  char *var;
  int ivar;

  double deadMass;

  void death();
};

}

#endif
#endif
