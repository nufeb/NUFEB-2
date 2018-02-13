/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/growth/energy,FixKineticsEnergy)

#else

#ifndef SRC_FIX_KINETICSENERGY_H
#define SRC_FIX_KINETICSENERGY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsEnergy : public Fix {
 public:
  FixKineticsEnergy(class LAMMPS *, int, char **);
  ~FixKineticsEnergy();
  void init();
  int setmask();
  void growth(double, int);

 private:
  char **var;
  int *ivar;

  double *radius;
  double *rmass;
  double *outerMass;
  double *outerRadius;

  int *mask;
  int *type;
  int nlocal;

  int nnus;                         // # of nutrients
  int ntypes;                       // # of species
  int epsflag;

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // computational domain size
  int nx, ny, nz;
  double vol;                       // grid volume and gas volume
  double temp;                      // gas transfer constant and temperature
  double EPSdens;                   // EPS density
  double *maintain;
  double *decay;
  double **DGRCat;

  double **catCoeff;                 // catabolism coefficients of species
  double **anabCoeff;                // anabolism  coefficients of species
  double **gYield;                   // yield coefficients
  double **gMonod;
  bool *nuConv;
  //double **minCatMonod;

  double **nuS;                    // nutrient concentration for all grids
  double **nuR;                    // nutrient consumption for all grids

  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class BIO *bio;

 // double minimal_monod(int, int, int);
  double grid_monod(int, int, int);
  void bio_update(double, int);
  double  biomass(int);
};

}

#endif
#endif

