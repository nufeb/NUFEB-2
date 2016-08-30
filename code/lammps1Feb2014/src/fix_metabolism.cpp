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

#include "fix_metabolism.h"
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
#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/KroneckerProduct>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;
using namespace Eigen;


/* ---------------------------------------------------------------------- */

FixMetabolism::FixMetabolism(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if (narg != 10) error->all(FLERR,"Not enough arguments in fix metabolism command");

  nevery = force->inumeric(FLERR,arg[3]);
  diffevery = force->inumeric(FLERR,arg[4]);
  if (nevery < 0 || diffevery < 0) error->all(FLERR,"Illegal fix growth command");

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[5+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[5+i][2]);
  }

  nx = atoi(arg[6]);
  ny = atoi(arg[7]);
  nz = atoi(arg[8]);

  if(strcmp(arg[9], "dirich") == 0) bflag = 1;
  else if(strcmp(arg[9], "neu") == 0) bflag = 2;
  else if(strcmp(arg[9], "mixed") == 0) bflag = 3;
  else error->all(FLERR,"Illegal boundary condition command");

//  // test
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

FixMetabolism::~FixMetabolism()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixMetabolism::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetabolism::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix  does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix  is invalid style");
  }

  diffT = input->variable->compute_equal(ivar[0]);
  nNus = atom->nNutrients;
  subname = atom->nuName;
  nuConc = atom->nuConc;
  catCoeff = atom->catCoeff;
  anabCoeff = atom->anabCoeff;
  diffCoeff = atom->diffCoeff;

  laplacian_matrix();
}

void FixMetabolism::laplacian_matrix()
{
  VectorXi ex1(nx);
  ex1.setOnes();
  MatrixXi  ex (nx, 3);
  ex << ex1, -3*ex1, ex1;
  VectorXi dx(3);
  dx << -1, 0, 1;
  SparseMatrix<int> Lx = spdiags(ex, dx, nx, ny, 3);

  VectorXi ey1(nx);
  ey1.setOnes();
  MatrixXi  ey (nx, 3);
  ey << ey1, -3*ey1, ey1;
  VectorXi dy(3);
  dy << -1, 0, 1;
  SparseMatrix<int> Ly = spdiags(ex, dy, nx, ny, 3);

  SparseMatrix<int> Ix(nx, nx);
  Ix.setIdentity();
  SparseMatrix<int> Iy(ny, ny);
  Iy.setIdentity();

  SparseMatrix<int> L2a = kroneckerProduct(Iy, Lx);
  SparseMatrix<int> L2b = kroneckerProduct(Ly, Ix);
  SparseMatrix<int> L2 = L2a + L2b;

  int N=nx*ny*nz;
  VectorXi ez1(N);
  ez1.setOnes();
  MatrixXi  ez (N, 2);
  ez << ez1, ez1;
  VectorXi dz(2);
  dz << -nx*ny, nx*ny;
  SparseMatrix<int> L = spdiags(ez, dz, N, N, 2);

  SparseMatrix<int> Iz(nz, nz);
  Iz.setIdentity();
  SparseMatrix<int> Aa = kroneckerProduct(Iz, L2);
  SparseMatrix<int> A = Aa + L;

  cout << "done!!!" <<endl;
}

SparseMatrix<int> FixMetabolism::spdiags(MatrixXi& B, VectorXi& d, int m, int n, int size_d)
{
  SparseMatrix<int> A(m,n);

  for (int k = 0; k < size_d; k++) {
    int i_min = max(0, (int)(-d(k)));
    int i_max = min(m - 1, (int)(n - d(k) - 1));
    int B_idx_start = m >= n ? d(k) : 0;

    for (int i = i_min; i <= i_max; i++) {
     A.insert(i, (int)(i+d(k))) = B(B_idx_start + i, k);
    }
  }
  return A;
}


/* ---------------------------------------------------------------------- */

void FixMetabolism::pre_force(int vflag)
{

}

/* ---------------------------------------------------------------------- */

void FixMetabolism::metabolism()
{

}

