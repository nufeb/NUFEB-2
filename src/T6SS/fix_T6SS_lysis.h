/* ----------------------------------------------------------------------
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

FixStyle(nufeb/T6SS/lysis, FixT6SSLysis)

#else

#ifndef LMP_FIX_T6SS_LYSIS_H
#define LMP_FIX_T6SS_LYSIS_H

#include "fix_growth.h"
namespace LAMMPS_NS {

class FixT6SSLysis : public FixGrowth {
public:
  FixT6SSLysis(class LAMMPS *, int, char **);
  ~FixT6SSLysis();
  int setmask();

  virtual void update_atoms();
  virtual void update_cells();
protected:
  int isub;
  double lysis_rate;
  double rbod_yield;
};

} // namespace LAMMPS_NS

#endif
#endif

    /* ERROR/WARNING messages:
     */
