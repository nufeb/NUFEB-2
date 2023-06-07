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

#ifdef FIX_CLASS

FixStyle(nufeb/plasmid/partition,FixPlasmidPartition)

#else

#ifndef LMP_FIX_PLASMID_PARTITION_H
#define LMP_FIX_PLASMID_PARTITION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPlasmidPartition : public Fix {
 public:
  FixPlasmidPartition(class LAMMPS *, int, char **);
  ~FixPlasmidPartition();

  void init();
  int setmask();
  void grow_arrays(int);
  void biology_nufeb();
  void compute();
  double memory_usage();

 private:
  void partition(int);
  int check_nucleoid(int, int, double);
  void relocate_plm_x(int, int);

  int nucleoid_flag;
  int fila_max;             // maximum # of filaments
  int *nfilas;              // # of filaments
  int ***fila;              // filament defined as straight line
  double **tfila;           // filament formation time
  double tmax_fila, v_fila; // maximum filament formation time, filament formation rate
  double diff_coef;         // diffusion coefficient for Browian motion
  double dt;                // timestep for plasmid movement
  int seed;
  int divflag;

  class RanPark *random;
  class AtomVecBacillus *avec;
  class FixPropertyPlasmid *fix_plm;
  class FixDivideBacillus *fix_div;
  class FixDivideBacillusMinicell *fix_div_mini;
};

}

#endif
#endif
