/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/monod,FixKineticsMonod)

#else

#ifndef SRC_FIX_KINETICSMONOD_H
#define SRC_FIX_KINETICSMONOD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsMonod : public Fix {
 public:
  FixKineticsMonod(class LAMMPS *, int, char **);
  ~FixKineticsMonod();
  int setmask();
  void init();
  void pre_force(int);

 private:
  char **var;
  int *ivar;

  double *radius;
  double *rmass;
  double *outerMass;
  double *outerRadius;

  int *mask;
  int* type;
  int nlocal;
  int nall;
  int diffevery;

  int nnus;                         // # of nutrients
  int ntypes;                       // # of species
  int nx, ny, nz;                   // number of grids in x y z axis
  int ngrids;                       //# of grids

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // computational domain size
  double vol;                       // grid volume and gas volume
  double temp;                      // gas transfer constant and temperature
  double EPSdens;                   // EPS density
  double *maintain;
  double *decay;
  double **DGRCat;

  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **gYield;                   // yield coefficients
  double **gMonod;
  double **minCatMonod;

  double **nuS;                    // nutrient concentration for all grids
  double **nuR;                    // nutrient consumption for all grids


  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class FixDiffusion *diffusion;
  class BIO *bio;

  void monod();
 // double minimal_monod(int, int, int);
  double grid_monod(int, int, int);
  void bio_update(double, int);
  double growth(int);
  int position(int);

};

}

#endif
#endif

