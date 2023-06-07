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

FixStyle(nufeb/diffusion_coeff,FixDiffusionCoeff)

#else

#ifndef LMP_FIX_DIFFUSION_COEFF_H
#define LMP_FIX_DIFFUSION_COEFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDiffusionCoeff : public Fix {
 public:
  FixDiffusionCoeff(class LAMMPS *, int, char **);
  virtual ~FixDiffusionCoeff() {}
  int setmask();
  virtual void init();
  virtual void post_physics_nufeb();
  virtual void compute();

 protected:
  int isub;
  int coeff_flag;
  double ratio;
  double vol;
  double const_coeff;

  class FixDiffusionReaction *fix_diffusion;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
