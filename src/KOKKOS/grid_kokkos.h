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

#include "grid.h"
#include "kokkos_type.h"

#ifndef LMP_GRID_KOKKOS_H
#define LMP_GRID_KOKKOS_H

namespace LAMMPS_NS {

class GridKokkos : public Grid {
 public:
  DAT::tdual_int_1d k_mask;
  DAT::tdual_float_1d k_bulk;
  DAT::tdual_float_2d k_conc;
  DAT::tdual_float_2d k_reac;
  DAT::tdual_float_2d k_diff_coeff;
  DAT::tdual_float_2d k_dens;
  DAT::tdual_int_2d k_boundary;
  DAT::tdual_float_3d k_growth;

  GridKokkos(class LAMMPS *);
  ~GridKokkos();

  void sync(const ExecutionSpace, unsigned int);
  void modified(const ExecutionSpace, unsigned int);
  void sync_overlapping_device(const ExecutionSpace, unsigned int);
  void sync_modify(ExecutionSpace, unsigned int, unsigned int);

 private:
  virtual class GridVec *new_gvec(const char *, int, int &);
};

}

#endif
