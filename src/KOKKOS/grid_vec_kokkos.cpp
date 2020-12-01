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

#include "grid_vec_kokkos.h"
#include "grid_kokkos.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

GridVecKokkos::GridVecKokkos(LAMMPS *lmp) : GridVec(lmp)
{
  kokkosable = 1;
}

void GridVecKokkos::setup()
{
  GridVec::setup();
  gridKK->modified(Host, GMASK_MASK);
  gridKK->sync(Device, GMASK_MASK);
}
