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
#include <Eigen/Eigen>

using namespace Eigen;

namespace LAMMPS_NS {

class FixMetabolism : public Fix {
 public:
  FixMetabolism(class LAMMPS *, int, char **);
  ~FixMetabolism();
  int setmask();
  void init();
  void pre_force(int);

 private:
  int nnus;                     // # of nutrients
  int ntypes;
  int nx, ny, nz;
  int ngrids;
  double stepx, stepy, stepz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double vol;
  double uniGas, temp;

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double *yield;

  double** metCoeff;               //metabolism coefficients of species
  int** matCons;                    //endergonic components of metabolic matrix
  VectorXd* vecR, vecConc;

  void metCoeff_calculus();
  void metabolism();
  double minimal_monod(int, int);

};

}

#endif
#endif

