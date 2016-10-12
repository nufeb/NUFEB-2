/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <fix_diffusion.h>
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>
#include <iomanip>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;
using namespace Eigen;


/* ---------------------------------------------------------------------- */

FixDiffusion::FixDiffusion(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if (narg != 11) error->all(FLERR,"Not enough arguments in fix diffusion command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix diffusion command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  //set grid size
  domain->nx = atoi(arg[5]);
  domain->ny = atoi(arg[6]);
  domain->nz = atoi(arg[7]);

  //set boundary condition flag:
  //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  if(strcmp(arg[8], "pp") == 0) xbcflag = 0;
  else if(strcmp(arg[8], "dd") == 0) xbcflag = 1;
  else if(strcmp(arg[8], "nd") == 0) xbcflag = 2;
  else if(strcmp(arg[8], "nn") == 0) xbcflag = 3;
  else if(strcmp(arg[8], "dn") == 0) xbcflag = 4;
  else error->all(FLERR,"Illegal x-axis boundary condition command");

  if(strcmp(arg[9], "pp") == 0) ybcflag = 0;
  else if(strcmp(arg[9], "dd") == 0) ybcflag = 1;
  else if(strcmp(arg[9], "nd") == 0) ybcflag = 2;
  else if(strcmp(arg[9], "nn") == 0) ybcflag = 3;
  else if(strcmp(arg[9], "dn") == 0) ybcflag = 4;
  else error->all(FLERR,"Illegal y-axis boundary condition command");

  if(strcmp(arg[10], "pp") == 0) zbcflag = 0;
  else if(strcmp(arg[10], "dd") == 0) zbcflag = 1;
  else if(strcmp(arg[10], "nd") == 0) zbcflag = 2;
  else if(strcmp(arg[10], "nn") == 0) zbcflag = 3;
  else if(strcmp(arg[10], "dn") == 0) zbcflag = 4;
  else error->all(FLERR,"Illegal z-axis boundary condition command");

//  // test
//  cout << xbcflag << ybcflag << zbcflag << endl;
//  for (int i = 1; i < atom->ntypes+1; i++){
//    printf("Atom info:\n");
//    printf("name = %s \n", atom->typeName[i]);
//    printf("type = %i, growth = %f, ks = %f, yield = %f \n",i, atom->growth[i], atom->ks[i], atom->yield[i]);
//    printf("Catabolism coeffs:\n");
//    for (int ii = 1; ii < atom->nNutrients +1; ii++) {
//      printf("%i ", atom->catCoeff[i][ii]);
//    }
//    printf("\n");
//    printf("Anabolism coeffs:\n");
//    for (int ii = 1; ii < atom->nNutrients +1; ii++) {
//      printf("%i ", atom->anabCoeff[i][ii]);
//    }
//    printf("\n");
//    printf("\n");
//  }
//
//  for (int i = 1; i < atom->nNutrients+1; i++){
//    printf("Nutrient info:\n");
//    printf("name = %s \n", atom->nuName[i]);
//    printf("type = %i, diff coeffs = %e, Scell = %e, Sbc=%e \n", i, atom->diffCoeff[i], atom->nuConc[i][0], atom->nuConc[i][1]);
//  }

}

/* ---------------------------------------------------------------------- */

FixDiffusion::~FixDiffusion()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixDiffusion::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDiffusion::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix diffusion does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix diffusion is invalid style");
  }

  diffT = input->variable->compute_equal(ivar[0]);
  nnus = atom->nNutrients;
  nuConc = atom->iniS;
  diffCoeff = atom->diffCoeff;

  nx = domain->nx;
  ny = domain->ny;
  nz = domain->nz;

  ngrids=nx*ny*nz;

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  double gridx = (xhi-xlo)/nx;
  double gridy = (yhi-ylo)/ny;
  double gridz = (zhi-zlo)/nz;

  if (!isEuqal(gridx, gridy, gridz)) error->all(FLERR,"Grid is not cubic");
  grid = gridx;

  nuS = memory->grow(atom->nuS, nnus+1, ngrids, "atom:nuS");
  nuR = memory->grow(atom->nuR, nnus+1, ngrids, "atom:nuR");

  //inlet concentration
  for (int i = 1; i <= nnus; i++) {
    for (int j = 0; j < ngrids; j++) {
      nuS[i][j] = nuConc[i][0];
    }
  }

  LAP = laplacian_matrix();
}

/* ----------------------------------------------------------------------
  build laplacian matrix
------------------------------------------------------------------------- */

SparseMatrix<double> FixDiffusion::laplacian_matrix()
{
  VectorXi ex1(nx);
  ex1.setOnes();
  MatrixXi  ex (nx, 3);
  ex << ex1, -3*ex1, ex1;
  VectorXi dx(3);
  dx << -1, 0, 1;
  SparseMatrix<double> Lx = spdiags(ex, dx, nx, ny, 3);
  VectorXi ey1(nx);
  ey1.setOnes();
  MatrixXi  ey (nx, 3);
  ey << ey1, -3*ey1, ey1;
  VectorXi dy(3);
  dy << -1, 0, 1;
  SparseMatrix<double> Ly = spdiags(ex, dy, nx, ny, 3);
  SparseMatrix<double> Ix(nx, nx);
  Ix.setIdentity();
  SparseMatrix<double> Iy(ny, ny);
  Iy.setIdentity();

  SparseMatrix<double> L2a = kroneckerProduct(Iy, Lx);
  SparseMatrix<double> L2b = kroneckerProduct(Ly, Ix);
  SparseMatrix<double> L2 = L2a + L2b;

  VectorXi ez1(ngrids);
  ez1.setOnes();
  MatrixXi  ez (ngrids, 2);
  ez << ez1, ez1;
  VectorXi dz(2);
  dz << -nx*ny, nx*ny;
  SparseMatrix<double> L = spdiags(ez, dz, ngrids, ngrids, 2);

  SparseMatrix<double> Iz(nz, nz);
  Iz.setIdentity();
  SparseMatrix<double> Aa = kroneckerProduct(Iz, L2);

  SparseMatrix<double> A = Aa + L;

  return A;
}

/* ----------------------------------------------------------------------
  extract and create sparse band and diagonal matrices
------------------------------------------------------------------------- */

SparseMatrix<double> FixDiffusion::spdiags(MatrixXi& B, VectorXi& d, int m, int n, int size_d)
{
  SparseMatrix<double> A(m,n);

  for (int k = 0; k < size_d; k++) {
    int i_min = max(0, (int)(-d(k)));
    int i_max = min(m - 1, (int)(n - d(k) - 1));
    int B_idx_start = m >= n ? d(k) : 0;

    for (int i = i_min; i <= i_max; i++) {
     A.insert(i, (double)(i+d(k))) = B(B_idx_start + i, k);
    }
  }
  return A;
}

/* ---------------------------------------------------------------------- */

void FixDiffusion::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  diffusion();
}

/* ----------------------------------------------------------------------
  solve diffusion equations and metabolism
------------------------------------------------------------------------- */

void FixDiffusion::diffusion()
{
  bool convergence = false;
  double tol = 1e-4; // tolerance for convergence criteria for nutrient balance equation
  int iteration = 0;
  double* r = new double[nnus+1]();
  bool* isConv = new bool[nnus+1]();
  double* maxBC = new double[nnus+1]();
  VectorXd* vecS = new VectorXd[nnus+1];
  VectorXd* vecR = new VectorXd[nnus+1];

  nuS = atom->nuS;
  nuR = atom->nuR;

  //initialization
  for (int i = 1; i <= nnus; i++) {
    //convert concentration and consumption data types to Eigen vector type
    vecS[i] = Map<VectorXd> (nuS[i], ngrids, 1);
    vecR[i] = Map<VectorXd> (nuR[i], ngrids, 1);
    isConv [i] = false;
    //get maximum boundary condition conc value
    double bc[6];
    for (int j = 0; j < 6; j++ ) {
      bc[j] = nuConc[i][j+1];
    }
    maxBC[i] = *max_element(bc, bc+6);
    r[i] = diffCoeff[i]*diffT/(grid*grid);
  }

  VectorXd BC(ngrids);
  VectorXd RES(ngrids);

  while (!convergence) {
    iteration ++;

    for (int i = 1; i <= nnus; i++) {
      if (!isConv[i] && atom->nuType[i] == 0) {
        isConv[0] = false;
        xbcm = nuConc[i][1];
        xbcp = nuConc[i][2];
        ybcm = nuConc[i][3];
        ybcp = nuConc[i][4];
        zbcm = nuConc[i][5];
        zbcp = nuConc[i][6];
        BC = bc_vec(vecS[i], grid);
        RES = LAP * vecS[i] + BC + (diffT * vecR[i]);
        RES = r[i] * RES;
        vecS[i] = RES + vecS[i];
        for (size_t j = 0; j < vecS[i].size(); j++) {
          if (vecS[i][j] < 0) vecS[i][j] = 1e-16;
        }
        double max = RES.array().abs().maxCoeff();
        double ratio = max/maxBC[i];
        if (ratio < tol) isConv[i] = true;
      }
      else isConv[0] = true;
    }
    if (isConv[0]) break;
  }
  cout << "number of iteration: " << iteration << endl;

  //convert concentration vector into normal data type
  for (int i = 1; i <= nnus; i++) {
    Map<MatrixXd>(nuS[i], vecS[i].rows(), vecS[i].cols()) =  vecS[i];
  }

//    for (int i = 1; i <= nnus; i++) {
//      cout << atom->nuName[i] << endl;
//      if (atom->nuType[i] == 0){
//        for (int j = 0; j < ngrids; j++) {
//          cout << nuS[i][j] << " ";
//        }
//      }
//      cout << endl;
//    }

  delete [] isConv;
  delete [] r;
  delete [] maxBC;
  delete [] vecS;
  delete [] vecR;
}

/* ----------------------------------------------------------------------
  build matrix of boundary condition
------------------------------------------------------------------------- */

VectorXd FixDiffusion::bc_vec(VectorXd& S, double h)
{
  VectorXd B(ngrids);
  B.setZero();
  int i_m;
  int i_p;

  //X-AXIS SURFACE
  for(int i = 1; i < ny +1; i++) {
    for(int j = 1; j < nz +1; j++) {
      int k = 1+(i-1)*nx+(j-1)*nx*ny;

      i_m = k-1 ;
      i_p = k+nx-2;
      //cout << i_m << " ip " << h << endl;
      if (xbcflag == 0) {
        B(i_m) = B(i_m)+S(i_p); //BOTTOM, X_MINUS
        B(i_p) = B(i_p)+S(i_m); //TOP, X_PLUS
      }
      else if (xbcflag == 1) {
        B(i_m) = B(i_m)+2*xbcm-S(i_m); //BOTTOM, X_MINUS
        B(i_p) = B(i_p)+2*xbcp-S(i_p); //TOP, X_PLUS
      }
      else if (xbcflag == 2) {
        B(i_m) = B(i_m)-h*xbcm+S(i_m); //BOTTOM,X_MINUS
        B(i_p) = B(i_p)+2*xbcp-S(i_p); //TOP, X_PLUS
      }
      else if (xbcflag == 3) {
        B(i_m) = B(i_m)-h*xbcm+S(i_m); //BOTTOM,X_MINUS
        B(i_p) = B(i_p)+h*xbcp+S(i_p); //TOP,X_PLUS
      }
      else if (xbcflag == 4) {
        B(i_m) = B(i_m)+2*xbcm-S(i_m); //BOTTOM, X_MINUS
        B(i_p) = B(i_p)+h*xbcp+S(i_p); //TOP,X_PLUS
      }
    }
  }

  //Y-AXIS SURFACE
  for(int i = 1; i < nz +1; i++) {
    i_p = (i)*nx*ny-nx;
    for (i_m = (i-1)*nx*ny; i_m <= (i-1)*nx*ny+nx-1; i_m++) {
      //error ip im?
      if (ybcflag == 0) {
        B(i_m) = B(i_m)+ S(i_p); //BOTTOM, Y_MINUS
        B(i_p)=B(i_p)+S(i_m); //TOP, Y_PLUS
      }
      else if (ybcflag == 1) {
        B(i_m)=B(i_m)+2*ybcm-S(i_m); //BOTTOM, Y_MINUS
        B(i_p)=B(i_p)+2*ybcp-S(i_p); //TOP, Y_PLUS
      }
      else if (ybcflag == 2) {
        B(i_m)=B(i_m)-h*ybcm+S(i_m); //BOTTOM,Y_MINUS
        B(i_p)=B(i_p)+2*ybcp-S(i_p); //TOP, Y_PLUS
      }
      else if (ybcflag == 3) {
        B(i_m)=B(i_m)-h*ybcm+S(i_m); //BOTTOM,Y_MINUS
        B(i_p)=B(i_p)+h*ybcp+S(i_p); //TOP,Y_PLUS
      }
      else if (ybcflag == 4) {
        B(i_m)=B(i_m)+2*ybcm-S(i_m); //BOTTOM, Y_MINUS
        B(i_p)=B(i_p)+h*ybcp+S(i_p); //TOP,Y_PLUS
      }
      i_p++;
    }
  }

  //Z-AXIS SURFACE
  i_p = nx*ny*(nz-1);
  for (i_m = 0; i_m < nx*ny; i_m++) {

    //error ip im?
    if (zbcflag == 0) {
      B(i_m)=B(i_m)+S(i_p); //BOTTOM, Z_MINUS
      B(i_p)=B(i_p)+S(i_m); //TOP, Z_PLUS
    }
    else if (zbcflag == 1) {
      B(i_m)=B(i_m)+2*zbcm-S(i_m); //BOTTOM, Z_MINUS
      B(i_p)=B(i_p)+2*zbcp-S(i_p); //TOP, Z_PLUS
    }
    else if (zbcflag == 2) {
      B(i_m)=B(i_m)-h*zbcm+S(i_m); //BOTTOM,Z_MINUS
      B(i_p)=B(i_p)+2*zbcp-S(i_p); //TOP, Z_PLUS
    }
    else if (zbcflag == 3) {
      B(i_m)=B(i_m)-h*zbcm+S(i_m); //BOTTOM,Z_MINUS
      B(i_p)=B(i_p)+h*zbcp+S(i_p); //TOP,Z_PLUS
    }
    else if (zbcflag == 4) {
      B(i_m)=B(i_m)+2*zbcm-S(i_m); //BOTTOM, Z_MINUS
      B(i_p)=B(i_p)+h*zbcp+S(i_p); //TOP,Z_PLUS
    }
    i_p++;
  }

  return B;
}

/* ----------------------------------------------------------------------
  compare double values for equality
------------------------------------------------------------------------- */

bool FixDiffusion::isEuqal(double a, double b, double c)
{
    double epsilon = 0.00001;
    if ((fabs(a - b) > epsilon)|| (fabs(a - b) > epsilon) || (fabs(a - b) > epsilon))
      return false;

    return true;
}

