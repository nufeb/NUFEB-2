/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/ph,FixKineticsPH)

#else

#ifndef SRC_FIX_KINETICSPH_H
#define SRC_FIX_KINETICSPH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsPH : public Fix {
 public:
  FixKineticsPH(class LAMMPS *, int, char **);
  ~FixKineticsPH();
  void init();
  int setmask();
  void solve_ph();

 private:
  int nnus;                     // # of nutrients
  double rth, temp;             // universal gas constant (thermodynamics) and temperature

  double **nuS;                 // nutrient concentration [nutrient][grid]
  double **kEq;                // equilibrium constants [nutrient][4]
  int **nuChr;                 // charge [nutrient][5charges]
  double ***activity;          // [nutrient][5 charges][grids]

  class FixKinetics *kinetics;
  class BIO *bio;

  void init_keq();
  void output_data();
};

}

#endif
#endif

