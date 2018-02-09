/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/growth/monod,FixKineticsMonod)

#else

#ifndef SRC_FIX_KINETICSMONOD_H
#define SRC_FIX_KINETICSMONOD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsMonod : public Fix {
 public:
  FixKineticsMonod(class LAMMPS *, int, char **);
  ~FixKineticsMonod();
  void init();
  int setmask();
  void growth(double, int);

 private:
  char **var;
  int *ivar;

  int isub, io2, inh4, ino2, ino3;   // nutrient index
  int idead, ieps;

  int *species;                      // species index 0 = unknow, 1 = het, 2 = aob, 3 = nob, 4 = eps, 5 = dead
  double **xtype;                    // type density [type][grids]
  double ***growrate;

  int nnus;                         // # of nutrients
  int ntypes;                       // # of species

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // computational domain size
  int nx, ny, nz;
  double vol;                       // grid volume and gas volume
  double EPSdens;                   // EPS density
  double etaHET;

  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class BIO *bio;

  void init_param();

};

}

#endif
#endif

