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

FixStyle(nufeb/plasmid/replicate,FixPlasmidReplication)

#else

#ifndef LMP_FIX_PLASMID_REPLICATION_H
#define LMP_FIX_PLASMID_REPLICATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPlasmidReplication : public Fix {
 public:
  FixPlasmidReplication(class LAMMPS *, int, char **);
  ~FixPlasmidReplication();

  void init();
  int setmask();
  void grow_arrays(int);
  void biology_nufeb();
  void compute();
  double memory_usage();

 private:
  void replication(int);

  double **nproteins;    // number of initiator proteins
  double alpha;
  int mean_protein, init_protein;
  int seed;

  class RanPark *random;

  class FixPropertyPlasmid *fix_plm;
};

}

#endif
#endif
