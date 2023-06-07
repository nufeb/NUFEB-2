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

FixStyle(nufeb/plasmid/kill,FixPlasmidKill)

#else

#ifndef LMP_FIX_PLASMID_KILL_H
#define LMP_FIX_PLASMID_KILL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPlasmidKill : public Fix {
 public:
  FixPlasmidKill(class LAMMPS *, int, char **);
  ~FixPlasmidKill() {}
  int modify_param(int, char **);

  void init() {}
  int setmask();
  void biology_nufeb();
  void compute();

 private:
  int idead;
  class FixPropertyPlasmid *fix_plasmid;
};

}

#endif
#endif
