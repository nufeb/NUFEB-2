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

#ifndef LMP_GRID_MASK_H
#define LMP_GRID_MASK_H

#define EMPTY_MASK     0x00000000
#define ALL_MASK       0xffffffff

#define BOUNDARY_MASK  0x000000ff
#define X_NB_MASK      0x00000001 // X negative boundary
#define X_PB_MASK      0x00000002 // X positive boundary
#define Y_NB_MASK      0x00000004 // Y negative boundary
#define Y_PB_MASK      0x00000008 // Y positive boundary
#define Z_NB_MASK      0x00000010 // Z negative boundary
#define Z_PB_MASK      0x00000020 // Z positive boundary
#define CORNER_MASK    0x00000100
#define GHOST_MASK     0x00000200

#endif

