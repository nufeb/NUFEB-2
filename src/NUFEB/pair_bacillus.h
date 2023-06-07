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

#ifdef PAIR_CLASS

PairStyle(bacillus,PairBacillus)

#else

#ifndef LMP_PAIR_BACILLUS_H
#define LMP_PAIR_BACILLUS_H

#include "pair.h"
#include <cmath>
#include "atom_vec_bacillus.h"

namespace LAMMPS_NS {

class PairBacillus : public Pair {
 public:
  PairBacillus(class LAMMPS *);
  ~PairBacillus();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  virtual void kernel_force(double R, int itype, int jtype,
    double& energy, double& fpair);

  struct Contact {
    int i, j;  // body (i.e. atom) indices (not tags)
    double fx,fy,fz;   // unscaled cohesive forces at contact
    double xi[3];      // coordinates of the contact point on ibody
    double xj[3];      // coordinates of the contact point on jbody
    double separation; // contact surface separation
  };

 protected:
  double **k_n;       // normal repulsion strength
  double **k_na;      // normal attraction strength
  double c_n;         // normal damping coefficient
  double c_t;         // tangential damping coefficient
  double mu;          // normal friction coefficient during gross sliding
  double cutoff;      // cutoff for interaction
  double *maxrad;   // per-type maximum enclosing radius

  int nmax;

  class AtomVecBacillus *avec;

  void allocate();

  // sphere-sphere interaction
  void sphere_against_sphere(int ibody, int jbody, int itype, int jtype,
                             double delx, double dely, double delz, double rsq,
                             double** v, double** f, double radi,
			     double radj, int evflag);
  // sphere-rod interaction
  void sphere_against_rod(int ibody, int jbody, int itype, int jtype,
                           double** x, double** v, double** f, double** torque,
                           double** angmom, AtomVecBacillus::Bonus *&ibonus,
			   int evflag);

  // rod-rod interactions
  void rod_against_rod(int ibody, int jbody, int itype, int jtype,
		       double** x, double** v, double** f, double** torque,
		       double** angmom, AtomVecBacillus::Bonus *&ibonus,
		       AtomVecBacillus::Bonus *&jbonus, int &num_contacts,
		       Contact &contact_list, double &evdwl, double* facc);

  // an edge vs. an rod from another body
  void interaction_rod_to_rod(int ibody, int edge_index_i, double* xmi,
                               double rounded_radius_i, int jbody, int edge_index_j,
                               double* xmj, double rounded_radius_j,
                               int itype, int jtype, double cutoff);

  // compute contact forces if contact points are detected
  void contact_forces(int ibody, int jbody, double *xi, double *xj,
		      double delx, double dely, double delz, double fx,
		      double fy, double fz, double** x, double** v,
		      double** angmom, double** f, double** torque, double* facc);

  // compute force and torque between two bodies given a pair of interacting points
  void pair_force_and_torque(int ibody, int jbody, double* pi, double* pj,
                             double r, double contact_dist, int itype, int jtype,
                             double** x, double** v, double** f, double** torque,
                             double** angmom, int jflag, double& energy, double* facc);

  // rescale the cohesive forces if a contact area is detected
  void rescale_cohesive_forces(double** x, double** f, double** torque,
                               Contact &contact_list, int &num_contacts,
                               int itype, int jtype, double* facc);

  // accumulate torque to a body given a force at a given point
  void sum_torque(double* xm, double *x, double fx, double fy, double fz, double* torque);

  // shortest distance between two rods (line segments)
  void distance_bt_rods(const double* x1, const double* x2,
                      const double* x3, const double* x4,
                      double* h1, double* h2, double& t1, double& t2, double& r);

  // shortest distance between sphere (point) and rod (line segments)
  void distance_bt_pt_rod(const double* q,
       const double* xi1, const double* xi2, double* h, double& d, double& t);

  void total_velocity(double* p, double *xcm, double* vcm, double *angmom,
                      double *inertia, double *quat, double* vi);

  double dot(double x1, double x2, double x3, double y1, double y2, double y3) { return x1*y1 + x2*y2 + x3*y3; }
  double norm(double x1, double x2, double x3) { return sqrt(dot(x1, x2, x3, x1, x2, x3)); }
  double dist(double u1, double u2, double u3, double v1, double v2, double v3) { return norm(u1-v1, u2-v2, u3-v3); }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair body/rounded/polyhedron requires atom style body rounded/polyhedron

Self-explanatory.

E: Pair body requires body style rounded/polyhedron

This pair style is specific to the rounded/polyhedron body style.

*/
