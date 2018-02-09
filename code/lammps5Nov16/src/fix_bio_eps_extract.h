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

FixStyle(eps_extract,FixEPSExtract)

#else

#ifndef LMP_FIX_EPS_EXTRACT_H
#define LMP_FIX_EPS_EXTRACT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEPSExtract : public Fix {
 public:
 
  FixEPSExtract(class LAMMPS *, int, char **);
  ~FixEPSExtract();
  int setmask();
  void init();
  void post_integrate();

 private:

  double growthFactor;
  int seed;
  char **var;
  int *ivar;
  int typeEPS;
  tagint maxtag_all;
  double xlo,xhi,ylo,yhi,zlo,zhi;

  class RanPark *random;
  class AtomVecBio *avec;
  class FixFluid *nufebFoam;
};
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use dynamic group with fix adapt atom

This is not yet supported.

E: Variable name for fix adapt does not exist

Self-explanatory.

E: Variable for fix adapt is invalid style

Only equal-style variables can be used.

E: Fix adapt pair style does not exist

Self-explanatory

E: Fix adapt pair style param not supported

The pair style does not know about the parameter you specified.

E: Fix adapt type pair range is not valid for pair hybrid sub-style

Self-explanatory.

E: Fix adapt kspace style does not exist

Self-explanatory.

E: Fix adapt requires atom attribute diameter

The atom style being used does not specify an atom diameter.

E: Fix adapt requires atom attribute charge

The atom style being used does not specify an atom charge.

E: Could not find fix adapt storage fix ID

This should not happen unless you explicitly deleted
a secondary fix that fix adapt created internally.

*/
