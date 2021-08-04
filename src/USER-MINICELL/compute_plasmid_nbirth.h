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

ComputeStyle(nufeb/plasmid/nbirth,ComputePlasmidNBirth)

#else

#ifndef LMP_COMPUTE_PLASMID_NBIRTH_H
#define LMP_COMPUTE_PLASMID_NBIRTH_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePlasmidNBirth : public Compute {
 public:
  ComputePlasmidNBirth(class LAMMPS *, int, char **);
  virtual ~ComputePlasmidNBirth() {}
  virtual void init() {}
  virtual double compute_scalar();
  void set_arrays(int);

  class AtomVecBacillus *avec;

 protected:
  int ncell, nplm;
  class FixPropertyPlasmid *fix_plasmid;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
