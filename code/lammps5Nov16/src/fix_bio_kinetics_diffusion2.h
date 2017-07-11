/*
 * fix_kinetics/diffusion2.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/diffusion2,FixKineticsDiffusion2)

#else

#ifndef SRC_FIX_KINETICS_DIFFUSION2_H
#define SRC_FIX_KINETICS_DIFFUSION2_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsDiffusion2 : public Fix {

 public:
  FixKineticsDiffusion2(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion2();
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
  int rstep;                    // steps skipped between Si+n-Si
  int rflag;
  double **nuR;
  double **nuS;
  double *diffD;

  double bl;
  double q, rvol, af;
  bool reactor;
  double *nuBS;

  int nx, ny, nz;               // # of non-ghost grids in x, y and z
  int nX, nY, nZ;               // # of all grids in x, y and z
  int nXYZ;                   // total # of grids
  double diffT;
  double xlo,xhi,ylo,yhi,zlo,zhi,izhi;
  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; // inlet BC concentrations for each surface

  double **xGrid;
  double **nuGrid;
  bool *ghost;

  class FixKinetics *kinetics;
  class BIO *bio;
  class AtomVecBio *avec;

  void update_grids(int);
  void compute_bc(double &, double *, int, double);
  void compute_bulk();
  void compute_bl();
  void compute_flux(double, double &, double *, double, int);
  bool isEuqal(double, double, double);
  double getMaxHeight();
  void test();
};

}

#endif
#endif

