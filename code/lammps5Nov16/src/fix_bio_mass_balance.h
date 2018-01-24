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

FixStyle(mbalance,FixMassBalance)

#else

#ifndef LMP_FIX_MASSBALANCE_H
#define LMP_FIX_MASSBALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMassBalance : public Fix {
 public:
  FixMassBalance(class LAMMPS *, int, char **);
 ~FixMassBalance();
  void init();
  int setmask();
  void end_of_step();

 private:

  int nlocal;
  int nall;
  int nnus;                         // # of nutrients

  int nevery;
  int cflag, nflag, mflag;

  double* pre_bmass;                  // total biomass in previoius step
  double* bmass;                     // total biomass in current step
  double pre_co2_carbon;
  double co2_carbon;
  double pre_nh3mass;
  double nh3mass;
  double pre_glu_carbon;
  double glu_carbon;
  double pre_o2mass;
  double o2mass;

  double **nuS;                    // nutrient concentration for all grids
  double **catCoeff;                 // catabolism coefficients of species
  double **anabCoeff;                // anabolism  coefficients of species
  double **gYield;                   // yield coefficients
  double vol;

  class FixKinetics *kinetics;
  class BIO *bio;

  void c_element_check();
  void n_element_check();
  void mass_check();

};

}

#endif
#endif


