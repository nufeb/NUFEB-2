/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(metabolism,FixMetabolism)

#else

#ifndef SRC_FIX_METABOLISM_H
#define SRC_FIX_METABOLISM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMetabolism : public Fix {
 public:
  FixMetabolism(class LAMMPS *, int, char **);
  ~FixMetabolism();
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
  double stepx, stepy, stepz;   // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;  //simulaton box size
  double vol, gvol;                //grid volume and gas volume
  double gasTran, temp;            //gas transfer constant and temperature

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double *yield;                   // yield coefficients

  double **metCoeff;               //metabolism coefficients of species
  int **matConsume;                //endergonic components of metabolic matrix
  double **nuS;                    //nutrient concentration for all grids
  double **nuR;                    //nutrient consumption for all grids

  void metCoeff_calculus();
  void metabolism();
  double minimal_monod(int, int);

};

}

#endif
#endif

