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
  char **subname;
  int diffevery;
  int outputevery;
//  double *xCell;
//  double *yCell;
//  double *zCell;
//  double *cellVol;
//  bool *ghost;

  double **nuConc;

  int **catCoeff;
  int **anabCoeff;

  double *diffCoeff;
  double diffT;

  int numCells;
  int nx, ny, nz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int bflag; // 1 = dirichlet, 2 = neumann, 3 = mixed
  double xstep, ystep, zstep;

  void laplacian_matrix();
  void metabolism();
  SparseMatrix<int> spdiags(MatrixXi&, VectorXi&, int, int, int);
 // arma::imat spdiags(const arma::imat&, const arma::irowvec&, int, int);
};

}

#endif
#endif

