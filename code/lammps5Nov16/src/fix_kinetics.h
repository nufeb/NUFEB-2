/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics,FixKinetics)

#else

#ifndef SRC_FIX_KINETICS_H
#define SRC_FIX_KINETICS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKinetics : public Fix {
  friend class FixKineticsMonod;
  friend class FixKineticsThermo;
  friend class FixDiffusion;
  friend class FixKineticsPH;

 public:
  FixKinetics(class LAMMPS *, int, char **);
  ~FixKinetics();
  int setmask();
  void init();

 private:
  char **var;
  int *ivar;

  int nx, ny, nz;                  // number of grids in x y z axis
  int ngrids;                      // # of grids

  double **nuS;                    // nutrient concentration [nutrient][grid]
  double **nuR;                    // nutrient consumption [nutrient][grid]
  double **nuGas;                  // gas correction [nutrient][grid]
  double **gYield;                 // inverse yield [type][grid]
  double **activity;               // activities of chemical species [nutrient][5 charges]
  double gVol, gasTrans;           // gas volume and gas transfer constant
  double temp, rth;                // uiversal gas constant (thermodynamics) and temperature

  class AtomVecBio *avec;
  class BIO *bio;
};

}

#endif
#endif

