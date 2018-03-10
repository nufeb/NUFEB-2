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

#ifdef COMPUTE_CLASS

ComputeStyle(roughness,ComputeNufebRough)

#else

#ifndef LMP_COMPUTE_ROUGH_H
#define LMP_COMPUTE_ROUGH_H

#include "compute.h"
#include "reduce_grid.h"

namespace LAMMPS_NS {

class ComputeNufebRough : public Compute, public ReduceGrid<ComputeNufebRough> {
  friend ReduceGrid<ComputeNufebRough>;

 public:
  ComputeNufebRough(class LAMMPS *, int, char **);
  virtual ~ComputeNufebRough();
  void init();
  virtual double compute_scalar();

 private:
  int nx, ny, nxy;
  double stepx, stepy;
  double *maxh;
  Grid<double, 2> grid;
  Subgrid<double, 2> subgrid;

  Subgrid<double, 2> get_subgrid() const { return subgrid; }
  template <typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      *result++ = maxh[*it];
    }
    return result;
  }
  template <typename InputIterator0, typename InputIterator1, typename BinaryOperation>
  InputIterator1 unpack_cells_reduce(InputIterator0 first, InputIterator0 last, InputIterator1 input, BinaryOperation op) {
    for (InputIterator0 it = first; it != last; ++it) {
      maxh[*it] = op(maxh[*it], *input++);
    }
    return input;
  }
  int get_cell_data_size(int n) const { return nxy; }
  bool is_bottom_most() const;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Temperature compute degrees of freedom < 0

This should not happen if you are calculating the temperature
on a valid set of atoms.

*/
