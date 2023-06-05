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

#ifndef LMP_GRID_VEC_H
#define LMP_GRID_VEC_H

#include <cstdio>
#include "pointers.h"

namespace LAMMPS_NS {

class GridVec : protected Pointers {
 public:
  int kokkosable;        // 1 if atom style is KOKKOS-enabled

  int nargcopy;          // copy of command-line args for atom_style command
  char **argcopy;        // used when AtomVec is realloced (restart,replicate)

  int size_forward;      // # of data in forward comm
  int size_exchange;     // # of data in exchange comm

  GridVec(class LAMMPS *);
  virtual ~GridVec() {}
  void store_args(int, char **);
  virtual void process_args(int, char **);
  virtual void init();
  virtual void grow(int) = 0;
  virtual void setup();

  virtual int pack_comm(int, int *, double *) = 0;
  virtual void unpack_comm(int, int *, double *) = 0;
  virtual int pack_exchange(int, int *, double *) = 0;
  virtual void unpack_exchange(int, int *, double *) = 0;
  
  virtual void set(int, char **) = 0;

 protected:
  int nmax;              // local copy of grid->nmax
};

}

#endif
