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
#define SRC_FIX_KINETICSMONOD _H

#include <fix.h>

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

  int nnus;                     // # of nutrients
  int ntypes;                   // # of species
  int nx, ny, nz;               // number of grids in x y z axis
  int ngrids;                   //# of grids

  double stepx, stepy, stepz;   // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;  //simulaton box size
  double vol, gvol;                //grid volume and gas volume
  double rg, temp;            //gas transfer constant and temperature

  int **matConsume;                // endergonic components of metabolic matrix
  double **catCoeff;               // catabolism coefficients of species
  double **anabCoeff;              // anabolism  coefficients of species
  double **metCoeff;               // metabolism coefficients of species
  double *yield;                   // yield coefficients

  double **nuS;                    //nutrient concentration for all grids
  double **nuR;                    //nutrient consumption for all grids

  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class BIO *bio;

  void create_metaMatrix();
  void monod();
  double minimal_monod(int, int);

};

}

#endif
#endif

