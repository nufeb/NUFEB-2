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

FixStyle(nufeb/death/plasmid,FixDeathPlasmid)

#else

#ifndef LMP_FIX_DEATH_PLASMID_H
#define LMP_FIX_DEATH_PLASMID_H

#include "fix_death.h"

namespace LAMMPS_NS {

class FixDeathPlasmid : public FixDeath {
 public:
  int compute_flag;

  FixDeathPlasmid(class LAMMPS *, int, char **);
  virtual ~FixDeathPlasmid() {}
  virtual void init() {}
  virtual void compute();

 private:
  int idead;
  class FixPropertyPlasmid *fix_plasmid;
};

}

#endif
#endif
