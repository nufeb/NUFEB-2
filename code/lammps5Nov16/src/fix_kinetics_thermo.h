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
  int setmask();
  void init();
  void pre_force(int);

 private:
  int nnus;                        // # of nutrients
  int ntypes;                      // # of species
  int nx, ny, nz;                  // number of grids in x y z axis
  int ngrids;                      //# of grids

  double rth, temp;                //Universal gas constant (thermodynamics) and temperature

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **nuGCoeff;               // Gibbs free energy coefficient [nutrient][5charges]
  double **typeGCoeff;             // Gibbs free energy coefficient [type][5charges]
  double *diss;                    // Gibbs free energy of dissipation
  double **iyield;                 // dynamic yield coeff [type][grid]

  double **nuS;                    // nutrient concentration for all grids
  double **dG0;
  double *khV;                     // Henry's constant
  int *liq2Gas;                    // liquids convert to gas

  double **DRGCat;                 // Gibbs free energy of catabolism [type][grid]
  double **DRGAn;                  // Gibbs free energy of anabolism [type][grid]

  class FixKinetics *kinetics;
  class BIO *bio;

  void thermo();
  void init_dG0();
  void init_KhV();
  void output_data();
};

}

#endif
#endif

