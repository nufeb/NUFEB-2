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
#ifdef INTEGRATE_CLASS
// clang-format off
IntegrateStyle(test/comm_grid,TestCommGrid)
// clang-format on
#else

#ifndef LMP_TEST_COMM_GRID_H
#define LMP_TEST_COMM_GRID_H

#include "integrate.h"

namespace LAMMPS_NS {

class TestCommGrid : public Integrate {
 public:
  TestCommGrid(class LAMMPS *, int, char **);
  ~TestCommGrid() {}
  void init() override;
  void setup(int flag) override;
  void setup_minimal(int) override;
  void run(int) override;
  void force_clear() {};
  void cleanup() {};

 private:
  bool check();
};

}

#endif
#endif
