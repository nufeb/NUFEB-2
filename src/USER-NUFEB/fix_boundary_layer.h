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

FixStyle(nufeb/boundary_layer,FixBoundaryLayer)

#else

#ifndef LMP_FIX_BOUNDARY_LAYER_H
#define LMP_FIX_BOUNDARY_LAYER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBoundaryLayer : public Fix {
 public:
  int compute_flag;

  FixBoundaryLayer(class LAMMPS *, int, char **);
  virtual ~FixBoundaryLayer() {}
  int setmask();
  int modify_param(int, char **);
  virtual void post_integrate();
  virtual void compute();

 protected:
  int xyl_flag, xyh_flag, yzl_flag, yzh_flag, xzl_flag, xzh_flag;
  int nlayers;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
