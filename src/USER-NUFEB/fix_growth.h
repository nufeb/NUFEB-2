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

#ifndef LMP_FIX_GROWTH_H
#define LMP_FIX_GROWTH_H

#include "fix.h"
#include "atom_vec_bacillus.h"

namespace LAMMPS_NS {

class FixGrowth : public Fix {
 public:
  FixGrowth(class LAMMPS *, int, char **);
  virtual ~FixGrowth() {}
  int modify_param(int, char **);

  virtual void chemistry_nufeb();
  virtual void biology_nufeb();

  int setmask();
  virtual void init();
  virtual void reset_dt();
  virtual void update_atoms() = 0;
  virtual void update_cells() = 0;
  
 protected:
  double dt;

  void update_atoms_coccus();
  void update_atoms_bacillus(AtomVecBacillus *&avec);
};

}

#endif

/* ERROR/WARNING messages:
*/
