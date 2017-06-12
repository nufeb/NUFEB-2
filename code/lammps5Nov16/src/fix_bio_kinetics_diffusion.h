/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/diffusion,FixKineticsDiffusion)

#else

#ifndef SRC_FIX_KINETICS_DIFFUSION_H
#define SRC_FIX_KINETICS_DIFFUSION_H

#include "fix.h"
#include "eigen/Eigen/Eigen"

using namespace Eigen;

namespace LAMMPS_NS {

class FixKineticsDiffusion : public Fix {

 public:
  FixKineticsDiffusion(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion();
  int setmask();
  void init();
  bool* diffusion(bool*, int, double);

 private:
  char **var;
  int *ivar;
  int nnus;                     // # of nutrients
  double stepx, stepy, stepz;

  int nlocal;
  int nall;
  double *rmass;
  int ntypes;                       // # of species
  int *mask;
  int *type;

  int sflag;
  double **iniS;                // inlet concentrations of nutrients
  double *diffCoeff;            // diffusion coefficients of nutrients
  double *mw;                   // molecular weights of nutrients
  double tol;                   // tolerance for convergence criteria for nutrient balance equation
  int rstep;                    // steps leave between Si+n-Si
  int rflag;
  double **nuR;
  double **nuS;
  double *r;
  double* maxBC;

  int nx, ny, nz;               // # of grids in x, y and z
  int ngrids;                   // total # of grids
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  double grid;                               // grid height
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; //inlet BC concentrations for each surface
  SparseMatrix<double> LAP;       //laplacian matrix
  SparseMatrix<double> I;       //sparse identity matrix
  VectorXd *nRES;

  class FixKinetics *kinetics;
  class BIO *bio;
  class AtomVecBio *avec;

  SparseMatrix<double> laplacian_matrix_3d();
  SparseMatrix<double> spdiags(MatrixXi&, VectorXi&, int, int, int);
//  template <class numeric_t>
//  SparseMatrix<numeric_t> spdiags2(const Matrix<numeric_t,-1,-1> &,
//            const VectorXi &, const int, const int);
  VectorXd bc_vec(VectorXd&, double);
  bool isEuqal(double, double, double);

  void test();
};

}

#endif
#endif

