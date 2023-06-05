/* ----------------------------------------------------------------------
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

FixStyle(nufeb/adhesion,FixAdhesion)

#else

#ifndef LMP_FIX_ADHESION_H
#define LMP_FIX_ADHESION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdhesion : public Fix {
 public:
  class NeighList *list;

  FixAdhesion(class LAMMPS *, int, char **);
  virtual ~FixAdhesion();
  void allocate();
  void init();
  void init_list(int, class NeighList *);
  int modify_param(int, char **);
  int setmask();
  virtual void post_force(int);

 protected:
  double **ah; //Hammaker constant
  double smin; //minimum separation
  double smax; //maximum seperatioin
  int allocated;

  template <int DISP>
  void compute();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
