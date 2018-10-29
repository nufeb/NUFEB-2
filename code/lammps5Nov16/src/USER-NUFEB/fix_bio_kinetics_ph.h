/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/ph,FixKineticsPH)

#else

#ifndef SRC_FIX_KINETICSPH_H
#define SRC_FIX_KINETICSPH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixKineticsPH : public Fix {
 public:
  FixKineticsPH(class LAMMPS *, int, char **);
  ~FixKineticsPH();
  void init();
  int setmask();
  void solve_ph(int first, int last); // first - first cell index
                                      // last - one past the last index (think stl algorithms)
  void buffer_ph();

 private:
  class FixKinetics *kinetics;
  class BIO *bio;

  double buffer_flag;

  void init_keq();
  void output_data();
  void compute_activity();
  void init_keq();
};

}

#endif
#endif

