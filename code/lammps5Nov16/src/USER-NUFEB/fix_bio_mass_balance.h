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

FixStyle(mbalance,FixVerify)

#else

#ifndef LMP_FIX_MASSBALANCE_H
#define LMP_FIX_MASSBALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixVerify : public Fix {
 public:
  FixVerify(class LAMMPS *, int, char **);
 ~FixVerify();
  void init();
  int setmask();
  void end_of_step();

 private:

  int nlocal;
  int nall;
  int nnus;                         // # of nutrients

  int nevery;
  int cflag, nflag, mflag;

  double total_bmass, pre_total_bmass;
  double co2_carbon, pre_co2_carbon;
  double glu_carbon, pre_glu_carbon;
  double nh3_nitrogen, pre_nh3_nitrogen;
  double no2_nitrogen, pre_no2_nitrogen;
  double no3_nitrogen, pre_no3_nitrogen;

  double **nuS;                    // nutrient concentration for all grids
  double **catCoeff;                 // catabolism coefficients of species
  double **anabCoeff;                // anabolism  coefficients of species
  double **gYield;                   // yield coefficients
  double vol;

  class FixKinetics *kinetics;
  class BIO *bio;

  void c_element_check();
  void nitrogen_mass_balance();
  void mass_check();

};

}

#endif
#endif


