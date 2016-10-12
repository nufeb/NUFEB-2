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
  int ntypes;
  int nx, ny, nz;
  int ngrids;
  double stepx, stepy, stepz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double vol, gvol;
  double gasTran, temp;

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double *yield;

  double **metCoeff;               //metabolism coefficients of species
  int **matConsume;                    //endergonic components of metabolic matrix
  double **nuS;
  double **nuR;

  void metCoeff_calculus();
  void metabolism();
  double minimal_monod(int, int);

};

}

#endif
#endif

