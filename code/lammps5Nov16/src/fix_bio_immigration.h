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

FixStyle(immigration,FixImmigration)

#else

#ifndef LMP_FIX_IMMIGRATION_H
#define LMP_FIX_IMMIGRATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixImmigration : public Fix {
 public:
	FixImmigration(class LAMMPS *, int, char **);
 ~FixImmigration();
  int setmask();
  void init();
  void post_integrate();

 private:
  char **var;
  int *ivar;

  int seed;
  double divMass, density;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  tagint maxtag_all;
  int zflag;

  double* gamfit(double*, int, double*);
  double find_z(double, double, double);
  void immgration();
  void find_maxid();
  char* create_type_name(char*);
  void test(int, int);

  class RanPark *random;
  class AtomVecBio *avec;
  class BIO *bio;
  class FixKinetics *kinetics;
};

}

#endif
#endif


