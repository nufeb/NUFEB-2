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

#include "fix_growth_imperf_denit_nitric_oxide.h"

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


#define FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE

/* ---------------------------------------------------------------------- */

FixGrowthImperfDenitNitricOxide::FixGrowthImperfDenitNitricOxide(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  printf("Found %d params\n",narg);
  if (narg != 31)
    error->all(FLERR, "Illegal fix nufeb/growth/ImperfDenitNO command. Expected 31  parameters, found ");

  if (!grid->chemostat_flag)
    error->all(FLERR, "fix nufeb/growth/ImperfDenitNO requires grid_style nufeb/chemostat");

  iss = -1;
  io2 = -1;
  ino3 = -1;
  ino2 = -1;
  ino = -1;
  
  k_s1 = 0.0;
  k_s2 = 0.0;
  k_s3 = 0.0;
  
  k_oh1 = 0.0;
  k_oh2 = 0.0;
  k_oh3 = 0.0;

  k_no3 = 0.0;
  k_no2 = 0.0;

  k_13no = 0.0;

  eta_g2 = 0.0;
  eta_g3 = 0.0;

  //TODO add default value init to eta_Y in full denit fix
  eta_Y = 0.0;

  growth = 0.0;
  yield = 1.0;
  decay = 0.0;

  iss = grid->find(arg[3]);
  if (iss < 0)
    error->all(FLERR, "Fix GrowthImperfDenitNO can't find substrate named " +std::string(arg[3]));
  k_s1 = utils::numeric(FLERR,arg[4],true,lmp);
  k_s2 = utils::numeric(FLERR,arg[5],true,lmp);
  k_s3 = utils::numeric(FLERR,arg[6],true,lmp);


  #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
  printf("Reading a fix/growth/ImperfDenitNO\n");
  printf("\tSubstrate: %s\n ", arg[3]);
  printf("\t\tk_s1: %E\n", k_s1);
  printf("\t\tk_s2: %E\n", k_s2);
  printf("\t\tk_s3: %E\n", k_s3);
  #endif

  io2 = grid->find(arg[7]);
  if (io2 < 0)
    error->all(FLERR, "Fix ImperfDenitNO can't find substrate named " +std::string(arg[7]));
  k_oh1 = utils::numeric(FLERR,arg[8],true,lmp);
  k_oh2 = utils::numeric(FLERR,arg[9],true,lmp);
  k_oh3 = utils::numeric(FLERR,arg[10],true,lmp);

  #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
  printf("\tSubstrate: %s\n ", arg[7]);
  printf("\t\tk_oh1: %E\n", k_oh1);
  printf("\t\tk_oh2: %E\n", k_oh2);
  printf("\t\tk_oh3: %E\n", k_oh3);
  #endif

  ino3 = grid->find(arg[11]);
  if (ino3 < 0)
    error->all(FLERR, "Fix ImperfDeniNO can't find substrate named " +std::string(arg[11]));
  k_no3 = utils::numeric(FLERR,arg[12],true,lmp);
  
  #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
  printf("\tSubstrate: %s\n ", arg[11]);
  printf("\t\tk_no3: %E\n", k_no3);
  #endif

  ino2 = grid->find(arg[13]);
  if (ino2 < 0)
    error->all(FLERR, "Fix ImperfDenitNO can't find substrate named " +std::string(arg[13]));
  k_no2 = utils::numeric(FLERR,arg[14],true,lmp);
  
  #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
  printf("\tSubstrate: %s\n ", arg[13]);
  printf("\t\tk_no2: %E\n", k_no2);
  #endif

  ino = grid->find(arg[15]);
  if (ino < 0)
    error->all(FLERR, "Fix ImperfDenitNO can't find substrate named " + std::string(arg[15]));
  k_13no = utils::numeric(FLERR,arg[16],true,lmp);

  #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
  printf("\tSubstrate: %s\n ", arg[15]);
  printf("\t\tk_13no: %E\n", k_13no);
  #endif

  inh = grid->find(arg[17]);
  if (inh < 0)
    error->all(FLERR, "Fix ImperfDenitNO can't find substrate named " + std::string(arg[17]));
  inxb = utils::numeric(FLERR,arg[18],true,lmp);

  #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
  printf("\tSubstrate: %s\n ", arg[17]);
  printf("\t\tinxb: %E\n", inxb);
  #endif

  int iarg = 19;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "growth") == 0) {
      growth = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
      printf("\tGrowth: %E\n ", growth);
      #endif
      iarg += 2;
    } else if (strcmp(arg[iarg], "yield") == 0) {
      yield = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
      printf("\tYield: %E\n ", yield);
      #endif
      iarg += 2;
    } else if (strcmp(arg[iarg], "decay") == 0) {
      decay = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
      printf("\tDecay: %E\n ", decay);
      #endif
      iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g2") == 0) {
       eta_g2 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
       printf("\teta_g2: %E\n ", eta_g2);
       #endif
       iarg += 2;
    } else if(strcmp(arg[iarg], "eta_g3") == 0) {
       eta_g3 = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
       printf("\teta_g3: %E\n ", eta_g3);
       #endif
       iarg += 2;
      } else if(strcmp(arg[iarg], "eta_Y") == 0) {
       eta_Y = utils::numeric(FLERR,arg[iarg+1],true,lmp);
       #ifdef FIX_GROWTH_IMPERF_DENIT_NO_VERBOSE
       printf("\teta_Y: %E\n ", eta_Y);
       #endif  
       iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/ImperfDenitNO command. Did not recognize argument name. Expected either growth, yield, decay,eta_Y, or eta_g2, eta_g3 got " + std::string(arg[iarg]));
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthImperfDenitNitricOxide::update_cells()
{
  double **conc = grid->conc;
  double **reac = grid->reac;
  double **dens = grid->dens;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      //using the terminology from Hiatt and Grady 2008
      //R1: aerobic growth 
      //R2: anoxic growth, nitrate -> nitrite
      //R3: anoxic growth, nitrite -> nitric oxide
      //the variable 'growth' here refers to mu_het, but is left as 'growth' within the class
      double mu = growth;
     
      // reusing a lot of concentrations, so for readability assign concentration at i to local vars 
      // compiler should optimize away under reasonable conditions (02, 03)
      double SS = conc[iss][i];
      double SO = conc[io2][i];
      double SNO3 = conc[ino3][i];
      double SNO2 = conc[ino2][i];
      double SNO = conc[ino][i];

      //TODO these rates are also calculated in update_atoms
      //should DRY
      //should also only calculate onece per timestep
      double r1 = mu * SS/(k_s1+SS) * SO/(k_oh1+SO);
      double r2 = mu * eta_g2 * SS/(k_s2+SS) * SNO3/(k_no3+SNO3) * k_oh2/(k_oh2 + SO); 
      //#TODO check on KOH3 vs KOH typo in original paper
      double r3 = mu * eta_g3 * (SS/(k_s3+SS)) * (SNO2/(k_no2+SNO2)) * (k_oh3/(k_oh3+SO)) * (k_13no/(k_13no+SNO));

      //TODO A and B can be calculated once at instantiation
      double A = (1-yield*eta_Y)/(1.143*yield*eta_Y);
      double B = (1-yield*eta_Y)/(0.571*yield*eta_Y);
      reac[iss][i] -= (1/yield *r1 + 1/(yield*eta_Y)*(r2+r3) ) * dens[igroup][i];
      reac[io2][i] -= (1-yield)/yield * (r1) * dens[igroup][i];
      reac[ino3][i] -= A * r2 * dens[igroup][i];
      reac[ino2][i] -= (-A*r2+B*r3) * dens[igroup][i];
      reac[ino][i] -= (-B*r3) * dens[igroup][i];
      reac[inh][i] -= inxb*(r1+r2+r3)*dens[igroup][i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthImperfDenitNitricOxide::update_atoms()
{
  double **conc = grid->conc;

  for (int i = 0; i < grid->ncells; i++) {
      //using the terminology from Hiatt and Grady 2008
      //R1: aerobic growth 
      //R2: anoxic growth, nitrate -> nitrite NO3 to NO2
      //R3: anoxic growth, nitrite -> nitric oxide NO2 to NO
      //the variable 'growth' here refers to mu_het, but is left as 'growth' within the class
      double mu = growth;
    
      // reusing a lot of concentrations, so for readability assign concentration at i to local vars 
      // compiler should optimize away under reasonable conditions (02, 03)
      double SS = conc[iss][i];
      double SO = conc[io2][i];
      double SNO3 = conc[ino3][i];
      double SNO2 = conc[ino2][i];
      double SNO = conc[ino][i];

      double r1 = mu * SS/(k_s1+SS) * SO/(k_oh1+SO);
      double r2 = mu * eta_g2 * SS/(k_s2+SS) * SNO3/(k_no3+SNO3) * k_oh2/(k_oh2 + SO); 
      //#TODO check on KOH3 vs KOH typo in original paper
      double r3 = mu * eta_g3 * (SS/(k_s3+SS)) * (SNO2/(k_no2+SNO2)) * (k_oh3/(k_oh3+SO)) * (k_13no/(k_13no+SNO));


      grid->growth[igroup][i][0] = r1 + r2 + r3 - decay;

  }

  update_atoms_coccus();
}
