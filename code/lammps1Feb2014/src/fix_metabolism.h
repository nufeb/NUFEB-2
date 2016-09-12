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
  int nNus;
  int diffevery;
  int outputevery;

  double **nuConc;

  int **catCoeff;
  int **anabCoeff;

  double *diffCoeff;
  double diffT;

  int numCells;
  int nx, ny, nz;
  int N;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int xbcflag, ybcflag, zbcflag; // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  double grid;
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp;

  SparseMatrix<double> A;

  SparseMatrix<double> laplacian_matrix();
  void metabolism();
  SparseMatrix<double> spdiags(MatrixXi&, VectorXi&, int, int, int);
  VectorXd bc_matrix(VectorXd&, double);
  VectorXd rs_matrix();
  bool isEuqal(double, double, double);

};

}

#endif
#endif

