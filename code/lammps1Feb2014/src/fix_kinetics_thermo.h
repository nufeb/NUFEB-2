/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/thermo,FixKineticsThermo)

#else

#ifndef SRC_FIX_KINETICSTHERMO_H
#define SRC_FIX_KINETICSTHERMO_H

#include <fix.h>

namespace LAMMPS_NS {

class FixKineticsThermo : public Fix {
 public:
  FixKineticsThermo(class LAMMPS *, int, char **);
  ~FixKineticsThermo();
  int setmask();
  void init();
  void pre_force(int);

 private:
  char **var;
  int *ivar;

  int nnus;                     // # of nutrients
  int ntypes;                   // # of species
  int nx, ny, nz;               // number of grids in x y z axis
  int ngrids;                   //# of grids

  double rth, temp;                //Universal gas constant (thermodynamics) and temperature

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **metCoeff;               // metabolism coefficients of species
  double **nuGCoeff;
  double **typeGCoeff;
  double *yield;                   // yield coefficients
  double **typeG;                  // type energy for all grids
  double *diss;                    // Gibbs free energy of dissipation
  double **iyield;
  int *ngflag;
  int *tgflag;

  double **nuS;                    //nutrient concentration for all grids
  double **nuG;                    //nutrient energy for all grids
  double **dG0;

  class FixKinetics *kinetics;
  class BIO *bio;

  void thermo();
  void init_dG0();
};

}

#endif
#endif

