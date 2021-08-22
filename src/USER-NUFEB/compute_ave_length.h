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

ComputeStyle(nufeb/ave_length,ComputeAveLength)

#else

#ifndef LMP_COMPUTE_AVE_LENGTH_H
#define LMP_COMPUTE_AVE_LENGTH_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAveLength : public Compute {
 public:
  ComputeAveLength(class LAMMPS *, int, char **);
  ~ComputeAveLength();
  void init() {}
  void compute_vector();

  class AtomVecBacillus *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
