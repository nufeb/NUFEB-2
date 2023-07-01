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

FixStyle(nufeb/ph,FixPH)

#else

#ifndef LMP_FIX_PH_H
#define LMP_FIX_PH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPH : public Fix {
 public:
  FixPH(class LAMMPS *, int, char **);
  virtual ~FixPH();

  int setmask();
  void chemistry_nufeb();
  void reactor_nufeb();
  virtual void init();
  virtual void compute();

  double *ph, *sh;           // per-grid attributes: pH and H+ concentration

 protected:
  double **sstc_gibbs;       // Gibbs energy for the substrate in 5 forms
                             // 0:Non-Hydrated, 1:Hydrated, 2:1st Deprotonation
                             // 3:2nd Deprotonation, 4:3rd Deprotonation
  int **ncharges;            // charge numbers for substrate in 5 forms
  int *form_id;              // substrate form for microbial utilisation

  double iph;                // initial ph
  double temp;               // temperature (K)

  int buff_flag;
  int ih, ina, icl;
  double phlo, phhi;

  double **keq;              // hydration and acid-base equilibrium constants
                             // 0: kd,  1:ka1,  2:ka2,  3:ka3

  double **act;             // activities of chemical forms for microbe uptaking
  double ***act_all;        // activities of 5 substrate forms

  void init_keq();
  void compute_ph(int, int);
  void compute_activity(int, int, double);
  void buffer_ph();
  void to_mol();
  void to_kg();
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
