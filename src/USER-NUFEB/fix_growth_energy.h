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

FixStyle(nufeb/growth/energy,FixGrowthEnergy)

#else

#ifndef LMP_FIX_GROWTH_ENERGY_H
#define LMP_FIX_GROWTH_ENERGY_H

#include "fix_growth.h"

namespace LAMMPS_NS {

class FixGrowthEnergy: public FixGrowth {
 public:
  FixGrowthEnergy(class LAMMPS *, int, char **);
  virtual ~FixGrowthEnergy();

  void biology_nufeb();
  void chemistry_nufeb();
  virtual void update_atoms();
  virtual void update_cells();

  virtual void init();

 protected:
  double *sub_gibbs;             // substrate Gibbs energy
  double *ks_coeff;              // ks coeffs
  double *cata_coeff;            // catabolic coeffs
  double *anab_coeff;            // anabolic coeffs
  double *decay_coeff;           // decay coeffs

  double uptake;                 // substrate update rate
  double decay;	                 // decay rate
  double maintain;               // maintenance rate
  double dissipation;            // dissipation energy
  double biomass_gibbs;          // biomass Gibbs energy
  double temp, gas_const;        // temperature (K) and ideal gas constant (KJ/mol.K)
  double mw_biomass;             // biomass molecular weight (g/mol)
  int e_donor;                    // electron donor

  double dgo_cata, dgo_anab;     // standard catabolic and anabolic Gibbs energy
  double *gibbs_cata;            // Gibbs free energy for catabolic reaction
  double *gibbs_anab;            // Gibbs free energy for anabolic reaction
  double *yield;
  double alfa, beta;

  class AtomVecBacillus *avec;

  void compute_dgo();
  void compute_dgr();
  double compute_monod(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:
*/
