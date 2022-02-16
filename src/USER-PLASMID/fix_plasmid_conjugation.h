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

FixStyle(nufeb/plasmid/conjugate,FixPlasmidConjugation)

#else

#ifndef LMP_FIX_PLASMID_CONJUGATION_H
#define LMP_FIX_PLASMID_CONJUGATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPlasmidConjugation : public Fix {
 public:
  FixPlasmidConjugation(class LAMMPS *, int, char **);
  ~FixPlasmidConjugation();

  void init();
  int setmask();
  void grow_arrays(int);
  void biology_nufeb();
  void compute();

 private:
  int irecipient, itrans;
  int seed;

  class RanPark *random;

  class FixPropertyPlasmid *fix_plm;
};

}

#endif
#endif
