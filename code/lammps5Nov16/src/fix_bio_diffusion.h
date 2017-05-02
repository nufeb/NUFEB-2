/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(diffusion,FixDiffusion)

#else

#ifndef SRC_FIX_DIFFUSION_H
#define SRC_FIX_DIFFUSION_H

#include "fix.h"
#include "Eigen/Eigen"

using namespace Eigen;

namespace LAMMPS_NS {

class FixDiffusion : public Fix {

 public:
  FixDiffusion(class LAMMPS *, int, char **);
  ~FixDiffusion();
  int setmask();
  void init();
  void diffusion();
  void consumption(VectorXd&, VectorXd&, int);

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
  double diffT;                 // diffusion timestamp
  double tol;                   // tolerance for convergence criteria for nutrient balance equation
  int rstep;                    // steps leave between Si+n-Si
  int rflag;
  double ** nuR;
  double ** nuS;
  double *r;
  double* maxBC;
  double *prevS;                 // maximum concentration in previous step
  double **gMonod;
  double **DGRCat;

  double *maintain;
  double *decay;
  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **gYield;                   // yield coefficients

  int nx, ny, nz;               // # of grids in x, y and z
  int ngrids;                   // total # of grids
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  double grid;                               // grid height
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; //inlet BC concentrations for each surface
  SparseMatrix<double> LAP;       //laplacian matrix
  SparseMatrix<double> I;       //sparse identity matrix

  class FixKinetics *kinetics;
  class BIO *bio;
  class AtomVecBio *avec;

  SparseMatrix<double> laplacian_matrix_3d();
  SparseMatrix<double> laplacian_matrix_2d();
  SparseMatrix<double> spdiags(MatrixXi&, VectorXi&, int, int, int);
  VectorXd bc_vec(VectorXd&, double);
  bool isEuqal(double, double, double);

  int position(int);
  double grid_monod(int, int, int, int, VectorXd&);

  void test();
};

}

#endif
#endif

