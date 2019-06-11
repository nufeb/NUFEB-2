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

FixStyle(nufeb/diffusion_reaction,FixDiffusionReaction)

#else

#ifndef LMP_FIX_DIFFUSION_REACTION_H
#define LMP_FIX_DIFFUSION_REACTION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDiffusionReaction : public Fix {
 public:
  bool compute_flag;

  FixDiffusionReaction(class LAMMPS *, int, char **);
  virtual ~FixDiffusionReaction();
  int setmask();
  int modify_param(int, char **);
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual double compute_scalar();
  virtual void reset_dt();
  
 protected:
  int isub;
  double diff_coef;
  int ncells;
  double *prev;
  double dt;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
