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

FixStyle(nufeb/eps_secretion,FixEPSSecretion)

#else

#ifndef LMP_FIX_EPS_SECRETION_H
#define LMP_FIX_EPS_SECRETION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEPSSecretion : public Fix {
 public:
  FixEPSSecretion(class LAMMPS *, int, char **);
  virtual ~FixEPSSecretion();

  int setmask();
  int modify_param(int, char **);
  virtual void biology_nufeb();
  virtual void post_neighbor();
  virtual void compute();
  
 protected:
  int type;
  int ieps;
  double ratio;
  double eps_density;
  int seed;

  class RanPark *random;
};
}

#endif
#endif

/* ERROR/WARNING messages:
*/
