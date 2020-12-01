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

#ifndef LMP_COMM_GRID_KOKKOS_H
#define LMP_COMM_GRID_KOKKOS_H

#include "comm_grid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class CommGridKokkos : public CommGrid {
 public:
  CommGridKokkos(class LAMMPS *);
  virtual ~CommGridKokkos();

  virtual void init();
  virtual void setup();                 // setup 3d comm pattern
  virtual void forward_comm();          // forward comm of grid data
  virtual void migrate();               // move cells to new procs

  template <class DeviceType>
  void forward_comm_device();
  
 protected:
  DAT::tdual_int_1d k_recv_begin, k_recv_end;
  DAT::tdual_int_1d k_send_begin, k_send_end;
  DAT::tdual_int_1d k_recv_cells, k_send_cells;
  DAT::tdual_xfloat_1d k_buf_recv, k_buf_send;

  DAT::tdual_int_1d k_recv_cells_self, k_send_cells_self;
  DAT::tdual_xfloat_1d k_buf_self;
  
  virtual void grow_recv(int);
  virtual void grow_send(int);
  virtual void grow_self(int);
};

}

#endif
