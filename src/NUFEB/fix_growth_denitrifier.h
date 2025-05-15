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

FixStyle(nufeb/growth/denit,FixGrowthDenit)

#else

#ifndef LMP_FIX_GROWTH_DENIT_H
#define LMP_FIX_GROWTH_DENIT_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthDenit: public FixGrowth {
 public:
  FixGrowthDenit(class LAMMPS *, int, char **);
  virtual ~FixGrowthDenit() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  int iss;
  int io2;
  int ino3;
  int ino2;
  int ino;
  int in2o;
  
  double k_s1;
  double k_s2;
  double k_s3;
  double k_s4;
  double k_s5;
  
  double k_oh1;
  double k_oh2;
  double k_oh3;
  double k_oh4;
  double k_oh5;

  double k_no3;
  double k_no2;
  double k_n2o;
  double k_no;

  double k_13no;
  double k_14no;
  double k_15no;

  double eta_g2;
  double eta_g3;
  double eta_g4;
  double eta_g5;

  double growth;
  double yield;
  double decay;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
