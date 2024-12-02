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

FixStyle(nufeb/T6SS/contact, FixT6SSContact)

#else

#ifndef LMP_FIX_T6SS_CONTACT_H
#define LMP_FIX_T6SS_CONTACT_H

#include "fix.h"
#include <unordered_map>
namespace LAMMPS_NS {

class FixT6SSContact : public Fix {
public:
  FixT6SSContact(class LAMMPS *, int, char **);
  ~FixT6SSContact();
  int setmask();

  virtual void biology_nufeb();
  
  void init_list(int, class NeighList *);
  void init();
  class NeighList *list;
private:
  void do_contact();
  
  int is_attack_possible(int i, int j);
  bool is_susceptible(int i, int j);
  int get_attack_index(int i);
  double periodic_displacement(int i, int j, int k);

  class RanPark *random;
  std::unordered_map<int, int> attacker_lookup;

  int ran_seed;
  int n_t6ss_types;
  int* t6ss_types;
  int* attacker_effectors;
  double* harpoon_lens;
  int* cooldowns;

  int n_susceptible_types;
  int* susceptible_types;
  int* susceptible_effectors;
  double* effector_success_probs;
  int* effector_groups;
  int* effector_effects;
  const int I_CANHIT_J = 1;
  const int J_CANHIT_I = 8;
  const int NO_ATTACKS = 0;
};

} // namespace LAMMPS_NS

#endif
#endif

    /* ERROR/WARNING messages:
     */
