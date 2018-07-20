/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/thermo,FixKineticsThermo)

#else

#ifndef SRC_FIX_KINETICSTHERMO_H
#define SRC_FIX_KINETICSTHERMO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsThermo : public Fix {
 public:
  FixKineticsThermo(class LAMMPS *, int, char **);
  ~FixKineticsThermo();
  void init();
  int setmask();
  void thermo(double dt);

 private:
  int yflag;                        // 0 = fixed yield 1 = dynamic yield
  int rflag;                       // 0 = open reactor 1 = closed reactor

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // simulaton box size
  double vol;                       // grid volume and gas volume

  double **dgzero;
  double *khv;                     // Henry's constant
  int *liqtogas;                    // liquids convert to gas

  class FixKinetics *kinetics;
  class BIO *bio;

  void init_dG0();
  void init_KhV();
};

}

#endif
#endif

