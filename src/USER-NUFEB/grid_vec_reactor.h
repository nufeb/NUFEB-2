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

#ifdef GRID_CLASS

GridStyle(nufeb/reactor,GridVecReactor)

#else

#ifndef LMP_GRID_VEC_REACTOR_H
#define LMP_GRID_VEC_REACTOR_H

#include "grid_vec.h"

namespace LAMMPS_NS {

class GridVecReactor : public GridVec {
 public:
  GridVecReactor(class LAMMPS *);
  ~GridVecReactor() {}
  void init();
  void grow(int);

  int pack_comm(int, int *, double *);
  void unpack_comm(int, int *, double *);
  int pack_exchange(int, int *, double *);
  void unpack_exchange(int, int *, double *);

  void set(int, char **);

 private:
  void set_reactor(int, double, double);

  int *mask;
  double **conc;    // concentration
  double **reac;    // reaction rate
  double **dens;    // density
  double ***growth; // growth rate
  double *bulk;     // bulk concentration
};

}

#endif
#endif
