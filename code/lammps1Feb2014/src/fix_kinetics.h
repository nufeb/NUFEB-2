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
  double temp;
  double **metCoeff;               //metabolism coefficients of species
  double **iyield;                  // inverse yield

  double **nuS;                    //nutrient concentration for all grids
  double **nuR;                    //nutrient consumption for all grids

  class AtomVecBio *avec;
  class BIO *bio;
};

}

#endif
#endif

