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
  char **var;
  int *ivar;
  int nNus;                     // # of nutrients
  int diffevery;
  int outputevery;

  double **nuConc;              // inlet concentrations of nutrients
  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double *yield;

  double *diffCoeff;            // diffusion coefficients of nutrients
  double diffT;                 // diffusion timestamp

  int nx, ny, nz;               // grid size
  int ngrids;                   // # of grids
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  double grid;                               // grid height
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; //inlet BC concentrations for each surface
  SparseMatrix<double> LAP;       //laplacian matrix

  SparseMatrix<double> laplacian_matrix();
  void diffusion();
  SparseMatrix<double> spdiags(MatrixXi&, VectorXi&, int, int, int);
  VectorXd bc_vec(VectorXd&, double);
  VectorXd _matrix();
  bool isEuqal(double, double, double);

  MatrixXd met_matrix();
  void removeRow(MatrixXd& , unsigned int);

};

}

#endif
#endif

