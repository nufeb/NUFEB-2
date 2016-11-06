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

#include <fix_kinetics.h>

namespace LAMMPS_NS {

class FixKineticsThermo : public FixKinetics {
 public:
  FixKineticsThermo(class LAMMPS *, int, char **);
  ~FixKineticsThermo();
  int setmask();
  void init();
  void pre_force(int);

 private:
  char **var;
  int *ivar;
  double temp;            //gas transfer constant and temperature

  void thermo();
};

}

#endif
#endif

