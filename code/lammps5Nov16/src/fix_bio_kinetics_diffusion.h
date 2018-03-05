/*
 * fix_kinetics/diffusionS.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/diffusion,FixKineticsDiffusion)

#else

#ifndef SRC_FIX_KINETICS_DIFFUSIONS_H
#define SRC_FIX_KINETICS_DIFFUSIONS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsDiffusion: public Fix {

 public:
  FixKineticsDiffusion(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion();

  char **var;
  int *ivar;
  int nvar;

  int nnus;                     // # of nutrients
  double stepx, stepy, stepz;

  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  int bulkflag;                           // 1=solve mass balance for bulk liquid

  double *rmass;
  int ntypes;                       // # of species
  int *mask;
  int *type;
  int shearflag, dragflag;

  double shearRate;
  double **iniS;                // inlet concentrations of nutrients
  double *diffCoeff;            // diffusion coefficients of nutrients
  double *mw;                   // molecular weights of nutrients
  double tol;                   // tolerance for convergence criteria for nutrient balance equation

  double **nuR;
  double **nuS;
  double *nuBS;                       // concentration in boundary layer [nutrient]
  double **nuGrid;                    // nutrient concentration in ghost mesh [nutrient][grid]
  double **xGrid;                     // grid coordinate [gird][3]
  bool *ghost;

  double *diffD;

  double bl;
  double q, rvol, af;
  bool reactor;
  int unit;                     // nutrient unit 0 = mol/l; 1 = kg/m3

  int nx, ny, nz;               // # of non-ghost grids in x, y and z
  int nX, nY, nZ;               // # of all grids in x, y and z
  int nXYZ;                   // total # of grids
  double diffT;
  double xlo, xhi, ylo, yhi, zlo, zhi, bzhi;
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; // inlet BC concentrations for each surface

  double **nuPrev;

  int recv_buff_size;
  int send_buff_size;
  double *recvbuff;
  double *sendbuff;
  int *convergences;

  MPI_Request *requests;
  MPI_Status *status;

  class FixKinetics *kinetics;
  class BIO *bio;
  class AtomVecBio *avec;

  int setmask();
  void init();
  int *diffusion(int*, int, double);
  void update_nuS();
  void update_grids();
  void compute_bc(double &, double *, int, double);
  void compute_bulk(int);
  void compute_bl();
  void compute_flux(double, double &, double *, double, int);
  bool isEuqal(double, double, double);
  int get_index(int);
};

}

#endif
#endif

