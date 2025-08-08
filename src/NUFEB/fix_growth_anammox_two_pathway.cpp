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

#include "fix_growth_anammox_two_pathway.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include "atom.h"
#include "error.h"
#include "grid.h"
#include "group.h"
#include "grid_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

//uncomment for more verbose output
//#define FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE

/* ---------------------------------------------------------------------- */

FixGrowthAnammoxTwoPathway::FixGrowthAnammoxTwoPathway(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg != 23)
    error->all(FLERR, "Illegal fix nufeb/growth/AnammoxTwoPathway command. Expected 23  parameters. ");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/AnammoxTwoPathway requires grid_style nufeb/chemostat");

  std::string name;
  int idx = 3;
  name = std::string(arg[idx]);
  io2 = grid->find(arg[idx++]);
  if (io2 < 0)
    error->all(FLERR, "Fix AnammoxTwoPathway can't find substrate named " + name);
  k_oh_an = utils::numeric(FLERR,arg[idx++],true,lmp);

  #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
  std::cout << "\tSubstrate: " << name  << std::endl;
  printf("\t\tk_oh_an: %E\n", k_oh_an);
  #endif

  //TODO DRY out parameter reading. Really should formalized and abstract out
  //the whole process at some point
  name = std::string(arg[idx]);
  ino2 = grid->find(arg[idx++]);
  if (ino2 < 0)
    error->all(FLERR, "Fix AnammoxTwoPathway can't find substrate named " + name);
  k_no2_an = utils::numeric(FLERR,arg[idx++],true,lmp);
  
  #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
  std::cout << "\tSubstrate: " << name  << std::endl;
  printf("\t\tk_no2_an: %E\n", k_no2_an);
  #endif

  name = std::string(arg[idx]);
  ino = grid->find(arg[idx++]);
  if (ino < 0)
    error->all(FLERR, "Fix AnammoxTwoPathway can't find substrate named " + name);
  k_no_an = utils::numeric(FLERR,arg[idx++],true,lmp);

  #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
  std::cout << "\tSubstrate: " << name  << std::endl;
  printf("\t\tk_no_an: %E\n", k_no_an);
  #endif

  name = std::string(arg[idx]);
  inh = grid->find(arg[idx++]);
  if (inh < 0)
    error->all(FLERR, "Fix AnammoxTwoPathway can't find substrate named " + name);
  k_nh_an = utils::numeric(FLERR,arg[idx++],true,lmp);
  inxb = utils::numeric(FLERR,arg[idx++],true,lmp);

  #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
  std::cout << "\tSubstrate: " << name << std::endl;
  printf("\t\tk_nh_an: %E\n", k_nh_an);
  printf("\t\tinxb: %E\n", inxb);
  #endif

  name = std::string(arg[idx]);
  ino3 = grid->find(arg[idx++]);
  if (ino3 < 0)
    error->all(FLERR, "Fix AnammoxTwoPathway can't find substrate named " + name);

  int iarg = idx;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      mu_max = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
      printf("\tGrowth: %E\n ", mu_max);
      #endif
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
      printf("\tYield: %E\n ", yield);
      #endif
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
      printf("\tDecay: %E\n ", decay);
      #endif
      iarg += 2;
    } else if(strcmp(arg[iarg], "eta_I_an") == 0) {
       eta_I_an = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
       printf("\teta_I_an: %E\n ", eta_I_an);
       #endif
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_S_an") == 0) {
       eta_S_an = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_ANAMMOX_TWO_PATHWAY_VERBOSE
       printf("\teta_S_an: %E\n ", eta_S_an);
       #endif
       iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/AnammoxTwoPathway command. Did not recognize argument name. Expected either growth, yield, decay, or eta_I_an, eta_S_an got " + std::string(arg[iarg]));
    }
  }
}

void FixGrowthAnammoxTwoPathway::computeRates(int cellIndex){
  double **conc = grid->conc;
  
  // reusing a lot of concentrations, so for readability assign concentration at i to local vars 
  // compiler should optimize away under reasonable conditions (02, 03)
  double SO = conc[io2][cellIndex];
  double SNO2 = conc[ino2][cellIndex];
  double SNO = conc[ino][cellIndex];
  double SNH = conc[inh][cellIndex];

  rI_AN = mu_max * eta_I_an * k_oh_an/(k_oh_an+SO) * SNO2/(k_no2_an+SNO2) * SNH/(k_nh_an+SNH);
  rS_AN = mu_max * eta_S_an * k_oh_an/(k_oh_an+SO) * SNO/(k_no_an+SNO) * SNH/(k_nh_an+SNH);
}

void FixGrowthAnammoxTwoPathway::update_cells()
{
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      computeRates(i);

      reac[ino3][i] +=  (rI_AN/1.14 + rS_AN/1.71) * dens[igroup][i];
      reac[ino2][i] -= (1/yield + 1/1.14) * rI_AN * dens[igroup][i];
      reac[ino][i] -= (1/yield + 1/1.71) * rS_AN * dens[igroup][i];
      reac[inh][i] -= ((1/yield-inxb)*rI_AN + (1/yield-inxb)*rS_AN) * dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthAnammoxTwoPathway::update_atoms()
{
  //TODO DRY with update_cells() and timestepping
  for (int i = 0; i < grid->ncells; i++) {
      computeRates(i);
      grid->growth[igroup][i][0] = rI_AN + rS_AN - decay;
  }

  update_atoms_coccus();
}
