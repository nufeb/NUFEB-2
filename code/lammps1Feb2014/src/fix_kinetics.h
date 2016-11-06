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
 public:
  FixKinetics(class LAMMPS *, int, char **);
  ~FixKinetics();
  int setmask();
  void init();

 private:
  char **var;
  int *ivar;

  int nnus;                     // # of nutrients
  int ntypes;                   // # of species
  int nx, ny, nz;               // number of grids in x y z axis
  int ngrids;                   //# of grids
  int xyz;
  double temp;

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **metCoeff;               //metabolism coefficients of species
  double *yield;                   // yield coefficients

  double **nuS;                    //nutrient concentration for all grids
  double **nuR;                    //nutrient consumption for all grids
  double **nuG;                    //energy for all grids

  class FixDiffusion *diffusion;
};

}

#endif
#endif

