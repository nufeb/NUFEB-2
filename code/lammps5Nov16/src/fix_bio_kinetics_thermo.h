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

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsThermo : public Fix {
 public:
  FixKineticsThermo(class LAMMPS *, int, char **);
  ~FixKineticsThermo();
  void init();
  int setmask();
  void thermo();

 private:
  char **var;
  int *ivar;

  int ntypes;                      // # of species
  int fixY;                        // 0 = fixed yield 1 = dynamic yield
  int closeR;                       // 0 = closed reactor 1 = open reactor

  double rth, temp;                //Universal gas constant (thermodynamics) and temperature
  double pressure;                 // Gas pressure

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // simulaton box size
  double vol;                       // grid volume and gas volume

  int nnus;                        // # of nutrients
  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **nuGCoeff;               // Gibbs free energy coefficient [nutrient][5charges]
  double **typeGCoeff;             // Gibbs free energy coefficient [type][5charges]
  double *diss;                    // Gibbs free energy of dissipation
  double *kLa;

  double **gYield;                 // dynamic yield coeff [type][grid]
  double **DRGCat;                 // Gibbs free energy of catabolism [type][grid]
  double **DRGAn;                  // Gibbs free energy of anabolism [type][grid]

  double **nuS;                    // nutrient concentration for all grids
  double **nuR;
  double ***activity;
  double **qGas;

  double **dG0;
  double *khV;                     // Henry's constant
  int *liq2Gas;                    // liquids convert to gas

  class FixKinetics *kinetics;
  class BIO *bio;

  void init_dG0();
  void init_KhV();
  //void output_data();
};

}

#endif
#endif

