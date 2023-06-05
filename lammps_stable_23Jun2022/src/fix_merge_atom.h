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

FixStyle(nufeb/merge_eps,FixMergeAtom)

#else

#ifndef LMP_FIX_MERGE_ATOM_H
#define LMP_FIX_MERGE_ATOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMergeAtom : public Fix {
 public:
  class NeighList *list;

  FixMergeAtom(class LAMMPS *, int, char **);
  virtual ~FixMergeAtom();
  void init();
  void init_list(int, class NeighList *);
  int setmask();
  virtual void biology_nufeb();
  virtual void post_neighbor();

 protected:
  double max_dia;
  double eps_den;
  int seed;

  void compute();
  class RanPark *random;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
