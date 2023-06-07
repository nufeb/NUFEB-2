/* ----------------------------------------------------------------------
   NUFEB - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nufeb/shear,FixShear)

#else

#ifndef LMP_FIX_SHEAR_H
#define LMP_FIX_SHEAR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixShear : public Fix {
 public:
  FixShear(class LAMMPS *, int, char **);
  virtual ~FixShear() {};
  int setmask();
  virtual void compute();
  virtual void post_force(int);

 protected:
  double shear_rate;
  double viscosity;
  double layer;

  int xflag, yflag;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/

