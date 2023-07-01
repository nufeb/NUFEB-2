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

FixStyle(nufeb/growth/cyano,FixGrowthCyano)

#else

#ifndef LMP_FIX_GROWTH_CYANO_H
#define LMP_FIX_GROWTH_CYANO_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthCyano: public FixGrowth {
 public:
  FixGrowthCyano(class LAMMPS *, int, char **);
  virtual ~FixGrowthCyano() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  int ilight;   // light
  int ico2;     // co2
  int igco2;    // gas co2
  int isuc;     // sucrose
  int io2;      // dissolved co2

  double light_affinity;
  double co2_affinity;

  double growth;
  double yield;
  double maintain;
  double decay;
  double suc_exp;
  
  class AtomVecBacillus *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
