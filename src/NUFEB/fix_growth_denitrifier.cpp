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

#include "fix_growth_denitrifier.h"

#include <cstdio>
#include <cstring>
#include <cmath>
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
//#define FIX_GROWTH_DENIT_VERBOSE

/* ---------------------------------------------------------------------- */

FixGrowthDenit::FixGrowthDenit(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg != 44)
    error->all(FLERR, "Illegal fix nufeb/growth/denit command. Expected 44  parameters. ");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/denit requires grid_style nufeb/chemostat");

  //TODO when we start refactoring parsing, look into making these static methods or some
  //other such thing which allows us to use const with parameters that don't change. also
  //sanity checking etc on input values.

  // a bunch of default value assignments were removed. don't put them back in
  // the philosophy is parsing should fail (maybe a bit more gracefully)
  // rather than feed in 'good' default parameters that silently make the sim
  // only *appear* to run the way the user expected

  iss = grid->find(arg[3]);
  if (iss < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[3]));
  k_s1 = utils::numeric(FLERR,arg[4],true,lmp);
  k_s2 = utils::numeric(FLERR,arg[5],true,lmp);
  k_s3 = utils::numeric(FLERR,arg[6],true,lmp);
  k_s4 = utils::numeric(FLERR,arg[7],true,lmp);
  k_s5 = utils::numeric(FLERR,arg[8],true,lmp);

  printf("Reading a fix/growth/denit\n");

  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("Reading a fix/growth/denit\n");
  printf("\tSubstrate: %s\n ", arg[3]);
  printf("\t\tk_s1: %E\n", k_s1);
  printf("\t\tk_s2: %E\n", k_s2);
  printf("\t\tk_s3: %E\n", k_s3);
  printf("\t\tk_s4: %E\n", k_s4);
  printf("\t\tk_s5: %E\n", k_s5);
  #endif

  io2 = grid->find(arg[9]);
  if (io2 < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[9]));
  k_oh1 = utils::numeric(FLERR,arg[10],true,lmp);
  k_oh2 = utils::numeric(FLERR,arg[11],true,lmp);
  k_oh3 = utils::numeric(FLERR,arg[12],true,lmp);
  k_oh4 = utils::numeric(FLERR,arg[13],true,lmp);
  k_oh5 = utils::numeric(FLERR,arg[14],true,lmp);

  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("\tSubstrate: %s\n ", arg[9]);
  printf("\t\tk_oh1: %E\n", k_oh1);
  printf("\t\tk_oh2: %E\n", k_oh2);
  printf("\t\tk_oh3: %E\n", k_oh3);
  printf("\t\tk_oh4: %E\n", k_oh4);
  printf("\t\tk_oh5: %E\n", k_oh5);
  #endif

  ino3 = grid->find(arg[15]);
  if (ino3 < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[15]));
  k_no3 = utils::numeric(FLERR,arg[16],true,lmp);
  
  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("\tSubstrate: %s\n ", arg[15]);
  printf("\t\tk_no3: %E\n", k_no3);
  #endif

  ino2 = grid->find(arg[17]);
  if (ino2 < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[17]));
  k_no2 = utils::numeric(FLERR,arg[18],true,lmp);
  
  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("\tSubstrate: %s\n ", arg[17]);
  printf("\t\tk_no2: %E\n", k_no2);
  #endif

  in2o = grid->find(arg[19]);
  if (in2o < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[19]));
  k_n2o = utils::numeric(FLERR,arg[20],true,lmp);

  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("\tSubstrate: %s\n ", arg[19]);
  printf("\t\tk_n2o: %E\n", k_n2o);
  #endif

  ino = grid->find(arg[21]);
  if (ino < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[21]));
  k_no = utils::numeric(FLERR,arg[22],true,lmp);
  k_I3no = utils::numeric(FLERR,arg[23],true,lmp);
  k_I4no = utils::numeric(FLERR,arg[24],true,lmp);
  k_I5no = utils::numeric(FLERR,arg[25],true,lmp);

  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("\tSubstrate: %s\n ", arg[21]);
  printf("\t\tk_no: %E\n", k_no);
  printf("\t\tk_I3no: %E\n", k_I3no);
  printf("\t\tk_I4no: %E\n", k_I4no);
  printf("\t\tk_I5no: %E\n", k_I5no);
  #endif

  inh = grid->find(arg[26]);
  if (inh < 0)
    error->all(FLERR, "Fix GrowthDenit can't find substrate named " + std::string(arg[26]));
  inxb = utils::numeric(FLERR,arg[27],true,lmp);

  #ifdef FIX_GROWTH_DENIT_VERBOSE
  printf("\tSubstrate: %s\n ", arg[26]);
  printf("\t\tinxb: %E\n", inxb);
  #endif

  int iarg = 28;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      mu_max = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_DENIT_VERBOSE
      printf("\tGrowth: %E\n ", mu_max);
      #endif
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_DENIT_VERBOSE
      printf("\tYield: %E\n ", yield);
      #endif
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_DENIT_VERBOSE
      printf("\tDecay: %E\n ", decay);
      #endif
      iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g2") == 0) {
       eta_g2 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_DENIT_VERBOSE
       printf("\teta_g2: %E\n ", eta_g2);
       #endif
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g3") == 0) {
       eta_g3 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_DENIT_VERBOSE
       printf("\teta_g3: %E\n ", eta_g3);
       #endif
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g4") == 0) {
       eta_g4 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_DENIT_VERBOSE
       printf("\teta_g4: %E\n ", eta_g4);
       #endif      
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g5") == 0) {
       eta_g5 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_DENIT_VERBOSE
       printf("\teta_g5: %E\n ", eta_g5);
       #endif  
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_Y") == 0) {
       eta_Y = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_DENIT_VERBOSE
       printf("\teta_Y: %E\n ", eta_Y);
       #endif  
       iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/denit command. Did not recognize argument name. Expected either growth, yield, decay,eta_Y, or eta_g2, through eta_g5");
    }
  }
  
  //these are used in reaction or yield results, but don't change between timesteps
  A = (1-yield*eta_Y)/(1.143*yield*eta_Y);
  B = (1-yield*eta_Y)/(0.571*yield*eta_Y);
}

//R1: aerobic growth 
double FixGrowthDenit::rate1(double SS, double SO)
{
      return (mu_max * SS/(k_s1+SS) * SO/(k_oh1+SO));
}

//R2: anoxic growth, nitrate -> nitrite
double FixGrowthDenit::rate2(double SS, double SNO3, double SO)
{
      return (mu_max * eta_g2 * SS/(k_s2+SS) * SNO3/(k_no3+SNO3) * k_oh2/(k_oh2 + SO));
}

//R3: anoxic growth, nitrite -> nitric oxide
double FixGrowthDenit::rate3(double SS, double SNO2, double SO, double SNO)
{
      return (mu_max * eta_g3 * (SS/(k_s3+SS)) * (SNO2/(k_no2+SNO2)) * (k_oh3/(k_oh3+SO)) * (k_I3no/(k_I3no+SNO)));
}

//R4: anoxic growth, nitric oxide -> nitrous oxide
double FixGrowthDenit::rate4(double SS, double SNO, double SO)
{
      return ( mu_max * eta_g4 * (SS/(k_s4+SS)) * (SNO/(k_no + SNO + (SNO*SNO)/k_I4no)) * (k_oh4/(k_oh4+SO)));
}

//R5: anoxic growth, nitrous oxide -> nitrogen
double FixGrowthDenit::rate5(double SS, double SN2O, double SO, double SNO)
{
      return (mu_max * eta_g5 * (SS/(k_s5+SS)) * (SN2O/(k_n2o + SN2O)) * (k_oh5/(k_oh5+SO)) * (k_I5no/(k_I5no+SNO)));
}

// extracted to DRY update_cells() and update_atoms()
// this may still get called twice on a timestep, but usually rarely and the
// add complexity of memoization doesn't seem worth it right now
void FixGrowthDenit::computeRates(int cellIndex){
  double **conc = grid->conc;
  
  // readability, compiler should optimize away under reasonable conditions (02, 03)
  double SS = conc[iss][cellIndex];
  double SO = conc[io2][cellIndex];
  double SNO3 = conc[ino3][cellIndex];
  double SNO2 = conc[ino2][cellIndex];
  double SNO = conc[ino][cellIndex];
  double SN2O = conc[in2o][cellIndex];

  r1 = rate1(SS, SO);
  r2 = rate2(SS, SNO3, SO);
  r3 = rate3(SS, SNO2, SO, SNO);
  r4 = rate4(SS, SNO, SO);
  r5 = rate5(SS, SN2O, SO, SNO);
}

void FixGrowthDenit::update_cells()
{
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++){
    if (grid->mask[i] & GRID_MASK){
      computeRates(i);

      reac[iss][i] -= ((1/yield * r1) + 1/(yield * eta_Y) * (r2+r3+r4+r5)) * dens[igroup][i];
      reac[io2][i] -= ((1-yield)/yield * r1) * dens[igroup][i];
      reac[ino3][i] -= A * r2 * dens[igroup][i];
      reac[ino2][i] += (A*r2 - B*r3) * dens[igroup][i];
      reac[ino][i] += (B*r3-B*r4) * dens[igroup][i];
      reac[in2o][i] += (B*r4-B*r5) * dens[igroup][i];
      reac[inh][i] -= inxb*(r1+r2+r3+r4+r5)*dens[igroup][i];
    }
  }
}

void FixGrowthDenit::update_atoms()
{
  //TODO DRY with respect to update_cells() and calc once per timestep
  for (int i = 0; i < grid->ncells; i++) {
      computeRates(i);
      grid->growth[igroup][i][0] = r1 + r2 + r3 + r4 +r5 - decay;
  }

  update_atoms_coccus();
}
