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

FixStyle(nufeb/death/diameter,FixDeathDiameter)

#else

#ifndef LMP_FIX_DEATH_DIAMETER_H
#define LMP_FIX_DEATH_DIAMETER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeathDiameter : public Fix {
 public:
  FixDeathDiameter(class LAMMPS *, int, char **);
  virtual ~FixDeathDiameter() {}
  int modify_param(int, char **);

  int setmask();
  virtual void init() {}
  virtual void biology_nufeb();
  virtual void compute();
  
 protected:
  int idead;
  int tdead;
  double diameter;
};

}

#endif
#endif
