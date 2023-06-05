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

#ifdef FIX_CLASS

FixStyle(nufeb/adhesion/bacillus,FixAdhesionBacillus)

#else

#ifndef LMP_FIX_ADHESION_BACILLUS_H
#define LMP_FIX_ADHESION_BACILLUS_H

#include "fix.h"
#include <cmath>
#include "atom_vec_bacillus.h"

namespace LAMMPS_NS {

class FixAdhesionBacillus : public Fix {
 public:
  class NeighList *list;

  FixAdhesionBacillus(class LAMMPS *, int, char **);
  virtual ~FixAdhesionBacillus() {}
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  virtual void post_force(int);
  void compute();
  // rod-rod interactions
  void rod_against_rod(int ibody, int jbody, double** x, double** v,
		       double** f, double** torque, double** angmom, double* rmass,
		       AtomVecBacillus::Bonus *&ibonus, AtomVecBacillus::Bonus *&jbonus);
  // shortest distance between two rods (line segments)
  void distance_bt_rods(const double* x1, const double* x2,
                      const double* x3, const double* x4,
                      double* h1, double* h2, double& t1, double& t2, double& r);
  // compute force and torque between two bodies given a pair of interacting points
  void adhesion_force_and_torque(int ibody, int jbody, double* pi, double* pj,
                             double r, double contact_dist, double** x, double** v,
			     double** f, double** torque, double** angmom, double* rmass);
  // accumulate torque to a body given a force at a given point
  void sum_torque(double* xm, double *x, double fx, double fy, double fz, double* torque);

  double dot(double x1, double x2, double x3, double y1, double y2, double y3) { return x1*y1 + x2*y2 + x3*y3; }
  double norm(double x1, double x2, double x3) { return sqrt(dot(x1, x2, x3, x1, x2, x3)); }

  double ke, cutoff;
  class AtomVecBacillus *avec;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
