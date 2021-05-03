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

FixStyle(nufeb/biomass_balance,FixBiomassBalance)

#else

#ifndef LMP_FIX_BIOMASS_BALANCE_H
#define LMP_FIX_BIOMASS_BALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBiomassBalance : public Fix {
 public:
  int compute_flag;

  FixBiomassBalance(class LAMMPS *, int, char **);
  virtual ~FixBiomassBalance();
  int setmask();
  int modify_param(int, char **);
  virtual void init();
  virtual void post_integrate();
  virtual void compute();

 protected:
  int nsubs;
  int *isub;
  double *inlet;

  double q;
  double rvol;
  double af;
  double xy;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
