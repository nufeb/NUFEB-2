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

FixStyle(nufeb/wall_adhesion,FixWallAdhesion)

#else

#ifndef LMP_FIX_WALL_ADHESION_H
#define LMP_FIX_WALL_ADHESION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallAdhesion : public Fix {
 public:
  FixWallAdhesion(class LAMMPS *, int, char **);
  virtual ~FixWallAdhesion() {}
  int setmask();
  virtual void post_force(int);
  virtual void compute();

 protected:
  double kn;
  int ieps;
  int eps_mask;
  int wallstyle;
  double lo,hi,cylradius;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
