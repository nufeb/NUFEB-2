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
#include "atom_vec_bacillus.h"

namespace LAMMPS_NS {

class FixPlasmidConjugation : public Fix {
 public:
  FixPlasmidConjugation(class LAMMPS *, int, char **);
  ~FixPlasmidConjugation() {};

  void init(){};
  int setmask();
  void biology_nufeb();
  void compute();
  void conjugate(int j, double* h1);
  // shortest distance between two rods (line segments)
  void distance_bt_rods(const double* x1, const double* x2,
                      const double* x3, const double* x4,
                      double* h1, double* h2, double& t1, double& t2, double& r);

 protected:
  int irecip, itrans;
  int irecipbit, itransbit;
  int ttrans, tflag;

  class FixPropertyPlasmid *fix_plm;
  class AtomVecBacillus *avec;

  void get_quat(double *, double *, double *);
  double dot(double x1, double x2, double x3, double y1, double y2, double y3) { return x1*y1 + x2*y2 + x3*y3; }
  double norm(double x1, double x2, double x3) { return sqrt(dot(x1, x2, x3, x1, x2, x3)); }
};

}

#endif
#endif
