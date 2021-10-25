/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nufeb/mass_transport,FixMassTransport)

#else

#ifndef LMP_FIX_MASS_TRANSPORT_H
#define LMP_FIX_MASS_TRANSPORT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMassTransport : public Fix {
 public:
  int compute_flag;

  FixMassTransport(class LAMMPS *, int, char **);
  ~FixMassTransport() {}
  int modify_param(int, char **);
  int setmask();
  void post_integrate();
  void compute();

 protected:
  int isub;
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
