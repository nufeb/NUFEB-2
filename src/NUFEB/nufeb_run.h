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

#ifdef INTEGRATE_CLASS

IntegrateStyle(nufeb,NufebRun)

#else

#ifndef LMP_NUFEB_RUN_H
#define LMP_NUFEB_RUN_H

#include "integrate.h"

namespace LAMMPS_NS {

class NufebRun : public Integrate {
 public:
  int ndiff;
  int npair;

  NufebRun(class LAMMPS *, int, char **);
  virtual ~NufebRun();
  virtual void init();
  virtual void setup(int flag);
  virtual void setup_minimal(int);
  virtual void run(int);
  virtual void reset_dt();
  void cleanup();

 protected:
  bool init_diff_flag, info;
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int torqueflag,extraflag;

  virtual void force_clear();

  double biodt;
  double diffdt;
  double difftol;
  int diffmax;
  double pairdt;
  double pairtol;
  int pairmax;
  int nfix_diffusion;

  class FixDensity *fix_density;
  class FixDiffusionReaction **fix_diffusion;
  class ComputePressure *comp_pressure;
  class ComputeKE *comp_ke;
  class ComputeVolume *comp_volume;

  FILE *profile;
  
  virtual void module_biology();
  virtual int module_chemistry();
  virtual double module_physics();
  virtual void module_post_physics();
  virtual void module_reactor();
  double get_time();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: No fixes defined, atoms won't move

If you are not using a fix like nve, nvt, npt then atom velocities and
coordinates will not be updated during timestepping.

E: KOKKOS package requires run_style verlet/kk

The KOKKOS package requires the Kokkos version of run_style verlet; the
regular version cannot be used.

*/
