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

#include "fix_T6SS_lysis.h"

#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "random_park.h"
#include "pair.h"
#include "update.h"
#include "grid.h"
#include "grid_masks.h"
#include "math_const.h"
#include <math.h>
#include <algorithm>
#include <string.h>
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixT6SSLysis::FixT6SSLysis(LAMMPS *lmp, int narg, char **arg)
    : FixGrowth(lmp, narg, arg){
    isub = -1;
    isub = grid->find(arg[3]);
    if(isub <0)
        error->all(FLERR,"Can't find substrate name");
    lysis_rate = utils::numeric(FLERR,arg[4],true,lmp);
    rbod_yield = utils::numeric(FLERR,arg[5],true,lmp);
 }

/* ---------------------------------------------------------------------- */
FixT6SSLysis::~FixT6SSLysis() { } 
/* ---------------------------------------------------------------------- */
void FixT6SSLysis::update_cells()
{
    double **conc = grid->conc;
    double **reac = grid->reac;
    double **dens = grid->dens;
    
    for (int i = 0;i < grid->ncells; i++){
        if(grid->mask[i] & GRID_MASK){
            reac[isub][i] += lysis_rate*dens[igroup][i]*rbod_yield;
        }
    }

}
void FixT6SSLysis::update_atoms()
{
    double **x = atom->x;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *biomass = atom->biomass;
    double *outer_radius = atom->outer_radius;
    double *outer_mass = atom->outer_mass;
    double **conc = grid->conc;
    
    double three_quarters_pi = (3.0 / (4.0 * MY_PI));
    const double four_thirds_pi = 4.0 * MY_PI / 3.0;
    const double third = 1.0 / 3.0;
    const double MIN_RAD = 1e-7;
    for(int i = 0; i < atom->nlocal; i++){
        if(atom->mask[i] & groupbit){
            const int cell = grid->cell(x[i]);
            const double density = rmass[i] /
                (four_thirds_pi * radius[i] * radius[i] * radius[i]);
            double min_mass = MIN_RAD * (four_thirds_pi * radius[i] * radius[i] * radius[i]);
            //printf("Bug %d radius %7g mass %7g density %7g min mass %7g\n",i,radius[i],rmass[i],density,min_mass);
            rmass[i] = rmass[i]*(1-lysis_rate*dt);
            if(rmass[i] < min_mass) rmass[i] = min_mass;
            radius[i] = pow(three_quarters_pi * (rmass[i]/density),third); 
            //printf("potentially new radius %7g mass %7g\n",radius[i],rmass[i]);
            if(radius[i] < MIN_RAD) radius[i] = MIN_RAD;
            //printf("new radius %7g mass %7g\n",radius[i],rmass[i]);
        }
    }
}
/* ---------------------------------------------------------------------- */

int FixT6SSLysis::setmask() {
  int mask = 0;
  mask |= CHEMISTRY_NUFEB;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

