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

FixStyle(nufeb/growth/anammox_two_pathway,FixGrowthAnammoxTwoPathway)

#else

#ifndef LMP_FIX_GROWTH_ANAMMOX_TWO_PATHWAY_H
#define LMP_FIX_GROWTH_ANAMMOX_TWO_PATHWAY_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthAnammoxTwoPathway: public FixGrowth {
 public:
  FixGrowthAnammoxTwoPathway(class LAMMPS *, int, char **);
  virtual ~FixGrowthAnammoxTwoPathway() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  int io2;
  int ino2;
  int ino;
  int inh;
  int ino3;
  
  double inxb;
  
  double k_oh_an;
  double k_no2_an;
  double k_nh_an;
  double k_no_an;
  
  double eta_I_an;
  double eta_S_an;

  double mu_max;
  double yield;
  double decay;

 private:
  void computeRates(int cellIndex);
  double rI_AN;
  double rS_AN;
};
}
#endif
#endif

/* ERROR/WARNING messages:
*/
