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

FixStyle(nufeb/growth/imperf_denit_nitric_oxide,FixGrowthImperfDenitNitricOxide)

#else

#ifndef LMP_FIX_GROWTH_IMPERF_DENIT_NO_H
#define LMP_FIX_GROWTH_IMPERF_DENIT_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthImperfDenitNitricOxide: public FixGrowth {
 public:
  FixGrowthImperfDenitNitricOxide(class LAMMPS *, int, char **);
  virtual ~FixGrowthImperfDenitNitricOxide() {}

  virtual void update_atoms();
  virtual void update_cells();

 protected:
  int iss;
  int io2;
  int ino3;
  int ino2;
  int ino;
  int in2o;
  int inh;
  
  double k_s1;
  double k_s2;
  double k_s3;
  double k_s4;
  double k_s5;
  
  double k_oh1;
  double k_oh2;
  double k_oh3;
  double k_oh4;
  double k_oh5;

  double k_no3;
  double k_no2;
  double k_n2o;
  double k_no;

  double k_13no;
  double k_14no;
  double k_15no;

  double inxb;

  double eta_g2;
  double eta_g3;
  double eta_g4;
  double eta_g5;

  double eta_Y;

  double mu_max;
  double yield;
  double decay;

private:
  // components of reaction or yield equations which can be reused and don't vary with timestep
  double A;
  double B;

  // rate equations
  //terminology from Hiatt and Grady 2008
  //R1: aerobic growth 
  //R2: anoxic growth, nitrate -> nitrite
  //R3: anoxic growth, nitrite -> nitric oxide
  //we use raten() to signify functions which calc RN
  //within the code we use rn=raten() - mainly for readability
  double rate1(double SS, double SO);
  double rate2(double SS, double SNO3, double SO);
  double rate3(double SS, double SNO2, double SO, double SNO);
  double rate4(double SS, double SNO, double SO);
  double rate5(double SS, double SN2O, double SO, double SNO);
 
  //rate for each of the above at a timestep and cell index
  //updated by computeRates
  //used in update_cells() and update_atoms()
  double r1;
  double r2;
  double r3;

  void computeRates(int cellIndex);
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
