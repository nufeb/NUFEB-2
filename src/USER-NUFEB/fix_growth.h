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

#ifndef LMP_FIX_MONOD_H
#define LMP_FIX_MONOD_H

#include "fix.h"
#include "atom_vec_bacillus.h"

namespace LAMMPS_NS {

class FixGrowth : public Fix {
 public:
  int compute_flag;
  int reaction_flag;
  int growth_flag;

  FixGrowth(class LAMMPS *, int, char **);
  virtual ~FixGrowth() {}
  int modify_param(int, char **);
  void update_atom();
  virtual void init();
  virtual void reset_dt();
  virtual int setmask();
  virtual void post_integrate();
  virtual void compute() = 0;
  
 protected:
  double dt;

  virtual void update_atoms() = 0;
  void update_atoms_coccus();
  void update_atoms_bacillus(AtomVecBacillus *&avec);
};

}

#endif

/* ERROR/WARNING messages:
*/
