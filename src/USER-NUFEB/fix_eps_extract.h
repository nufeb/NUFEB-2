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

FixStyle(nufeb/eps_extract,FixEPSExtract)

#else

#ifndef LMP_FIX_EPS_EXTRACT_H
#define LMP_FIX_EPS_EXTRACT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEPSExtract : public Fix {
 public:
  FixEPSExtract(class LAMMPS *, int, char **);
  ~FixEPSExtract();

  int setmask();
  int modify_param(int, char **);
  void biology_nufeb();
  void post_neighbor();
  void compute();
  
 private:
  int type;
  int ieps;
  double ratio;
  double density;
  int seed;

  class RanPark *random;
};
}

#endif
#endif

/* ERROR/WARNING messages:
*/
