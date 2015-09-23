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

FixStyle(diffnugrowth,FixDiffNuGrowth)

#else

#ifndef LMP_FIX_DIFFNUGROWTH_H
#define LMP_FIX_DIFFNUGROWTH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDiffNuGrowth : public Fix {
 public:
  FixDiffNuGrowth(class LAMMPS *, int, char **);
  ~FixDiffNuGrowth();
  int setmask();
  void init();
  void pre_force(int);

 private:

  char **var;
  int *ivar;
  double *xCell;
  double *yCell;
  double *zCell;
  double *cellVol;
  bool *boundary;
  int numCells;
  int nx, ny, nz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  bool xloBound, xhiBound, yloBound, yhiBound, zloBound, zhiBound;
  double xstep, ystep, zstep;
  void change_dia();
};

}

#endif
#endif
