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

#ifndef KOKKOS_BASE_H
#define KOKKOS_BASE_H

#include "kokkos_type.h"

namespace LAMMPS_NS {

class KokkosBase {
 public:
  KokkosBase() {}

  // Pair
  virtual int pack_forward_comm_kokkos(int, DAT::tdual_int_2d,
                                       int, DAT::tdual_xfloat_1d &,
                                       int, int *) {return 0;};
  virtual void unpack_forward_comm_kokkos(int, int, DAT::tdual_xfloat_1d &) {}

  // Region
  virtual void match_all_kokkos(int, DAT::tdual_int_1d) {}

  // Fix
  virtual int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
                                   DAT::tdual_int_1d k_sendlist,
                                   DAT::tdual_int_1d k_copylist,
                                   ExecutionSpace space, int dim,
                                   X_FLOAT lo, X_FLOAT hi) { return 0; }
  virtual void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                                      DAT::tdual_int_1d &indices,int nrecv,
                                      int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
                                      ExecutionSpace space) {}
};

}

#endif

/* ERROR/WARNING messages:

*/
