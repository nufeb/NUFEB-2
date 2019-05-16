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

#ifndef LMP_COMM_GRID_H
#define LMP_COMM_GRID_H

#include "pointers.h"

namespace LAMMPS_NS {

class CommGrid : protected Pointers {
 public:
  CommGrid(class LAMMPS *);
  virtual ~CommGrid();

  virtual void init();
  virtual void setup();                 // setup 3d comm pattern
  virtual void forward_comm();          // forward comm of grid data
  virtual void exchange();              // move cells to new procs
  
 protected:
  int size_forward;                     // # of data in forward comm
  int size_exchange;                    // # of data in exchange comm
  int max_size;                         // maximum size between forward and
                                        // exchange comm data

  int nrecv;
  int nsend;
  int nrecvproc;
  int nsendproc;
  int *recvproc;
  int *sendproc;
  int *recv_begin, *recv_end;
  int *send_begin, *send_end;
  int *recv_cells;
  int *send_cells;
  double *buf_recv;
  double *buf_send;

  int nrecv_self;
  int nsend_self;
  int *recv_cells_self;
  int *send_cells_self;
  double *buf_self;
  
  MPI_Request *requests;
  
  void grow_recv(int);
  void grow_send(int);
  void grow_self(int);
  int intersect(int *lo1, int *hi1,
		int *lo2, int *hi2,
		int ext1, int ext2,
		int xshift, int yshift, int zshift,
		int *intlo, int *inthi);
};

}

#endif
