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
  int isub;
  double diff_coeff;

  FixDiffusionReaction(class LAMMPS *, int, char **);
  virtual ~FixDiffusionReaction();
  int setmask() {return 0;}
  int modify_param(int, char **);
  void reset_dt();

  virtual void init();
  virtual double compute_scalar();
  virtual void compute_initial();
  virtual void compute_final();
  virtual void closed_system_initial();
  virtual void closed_system_scaleup(double);

 protected:
  int ncells;
  double *prev;		       // substrate concentration at n-1 step
  double dt;
  int *boundary;

  double *penult;	       // substrate concentration at n-2 step
  int closed_system;

};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
