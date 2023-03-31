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

#include <cstdio>
#include <cstring>
#include <cmath>
#include "fix_diffusion_reaction.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "grid.h"
#include "grid_masks.h"
#include "memory.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define THRESHOLD_CONC 1E-20

/* ---------------------------------------------------------------------- */

FixDiffusionReaction::FixDiffusionReaction(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR,"Illegal fix nufeb/diffusion_reaction command");

  dynamic_group_allow = 1;
  scalar_flag = 1;

  closed_system = 0;
  ncells = 0;
  dt = 1.0;
  prev = nullptr;
  penult = nullptr;
  boundary = nullptr;

  isub = grid->find(arg[3]);
  if (isub < 0)
    error->all(FLERR, "Can't find substrate for nufeb/diffusion_reaction");

  diff_coeff = utils::numeric(FLERR,arg[4],true,lmp);
}

/* ---------------------------------------------------------------------- */

FixDiffusionReaction::~FixDiffusionReaction()
{
  if (copymode) return;
  memory->destroy(prev);
  if (closed_system) memory->destroy(penult);
}

/* ---------------------------------------------------------------------- */

int FixDiffusionReaction::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "nevery") == 0) {
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery <= 0) error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::init()
{
  boundary = grid->boundary[isub];
  ncells = grid->ncells;
  prev = memory->create(prev, ncells, "nufeb/diffusion_reaction:prev");
  for (int i = 0; i < ncells; i++) {
    prev[i] = 0.0;
    grid->diff_coeff[isub][i] = diff_coeff;
  }

  int fix = 0;
  for (int i = 0; i < 6; i++) {
    if (boundary[i] == DIRICHLET) {
      fix = 1;
      break;
    }
  }
  if (!fix) closed_system = 1;

  if (closed_system) {
    penult = memory->create(penult, ncells, "nufeb/diffusion_reaction:penult");
    for (int i = 0; i < ncells; i++)
      penult[i] = 0.0;
  }
  dt = update->dt;

}

/* ---------------------------------------------------------------------- */

double FixDiffusionReaction::compute_scalar()
{
  double result = 0.0;
  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      double res = fabs((grid->conc[isub][i] - prev[i]) / prev[i]);
      if (closed_system) {
        double res2 = fabs((prev[i] - penult[i]) / penult[i]);
        res = fabs(res - res2);
      }
      result = MAX(result, res);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_MAX, world);
  return result;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::compute_initial()
{
  if (ncells != grid->ncells) {
    ncells = grid->ncells;
    prev = memory->grow(prev, ncells, "nufeb/diffusion_reaction:prev");
    if (closed_system) penult = memory->grow(penult, ncells, "nufeb/diffusion_reaction:penult");
  }

  for (int i = 0; i < grid->ncells; i++) {
    // Dirichlet boundary conditions
    if (grid->mask[i] & X_NB_MASK && (boundary[0] == DIRICHLET)) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & X_PB_MASK && boundary[1] == DIRICHLET) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Y_NB_MASK && boundary[2] == DIRICHLET) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Y_PB_MASK && boundary[3] == DIRICHLET) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Z_NB_MASK && boundary[4] == DIRICHLET) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & Z_PB_MASK && boundary[5] == DIRICHLET) {
      grid->conc[isub][i] = grid->bulk[isub];
    } else if (grid->mask[i] & BLAYER_MASK) {
      grid->conc[isub][i] = grid->bulk[isub];
    }
    grid->reac[isub][i] = 0.0;
  }
}



/* ---------------------------------------------------------------------- */

void FixDiffusionReaction::compute_final()
{
  int nx = grid->subbox[0];
  int nxy = grid->subbox[0] * grid->subbox[1];
  for (int i = 0; i < grid->ncells; i++) {
    // Neumann boundary conditions
    if (grid->mask[i] & X_NB_MASK && boundary[0] == NEUMANN) {
      grid->conc[isub][i] = grid->conc[isub][i+1];
    } else if (grid->mask[i] & X_PB_MASK && boundary[1] == NEUMANN) {
      grid->conc[isub][i] = grid->conc[isub][i-1];
    } else if (grid->mask[i] & Y_NB_MASK && boundary[2] == NEUMANN) {
      int py = i + nx;
      grid->conc[isub][i] = grid->conc[isub][py];
    } else if (grid->mask[i] & Y_PB_MASK && boundary[3] == NEUMANN) {
      int py = i - nx;
      grid->conc[isub][i] = grid->conc[isub][py];
    } else if (grid->mask[i] & Z_NB_MASK && boundary[4] == NEUMANN) {
      int pz = i + nxy;
      grid->conc[isub][i] = grid->conc[isub][pz];
    } else if (grid->mask[i] & Z_PB_MASK && boundary[5] == NEUMANN) {
      int pz = i - nxy;
      grid->conc[isub][i] = grid->conc[isub][pz];
    }
    if (closed_system) penult[i] = prev[i];
    prev[i] = grid->conc[isub][i];
  }

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      int nx = i - 1;
      int px = i + 1;
      int ny = i - grid->subbox[0];
      int py = i + grid->subbox[0];
      int nz = i - nxy;
      int pz = i + nxy;
      double dnx = grid->diff_coeff[isub][i] * (prev[i] - prev[nx]) / grid->cell_size;
      double dpx = grid->diff_coeff[isub][i] * (prev[px] - prev[i]) / grid->cell_size;
      double ddx = (dpx - dnx) / grid->cell_size;
      double dny = grid->diff_coeff[isub][i] * (prev[i] - prev[ny]) / grid->cell_size;
      double dpy = grid->diff_coeff[isub][i] * (prev[py] - prev[i]) / grid->cell_size;
      double ddy = (dpy - dny) / grid->cell_size;
      double dnz = grid->diff_coeff[isub][i] * (prev[i] - prev[nz]) / grid->cell_size;
      double dpz = grid->diff_coeff[isub][i] * (prev[pz] - prev[i]) / grid->cell_size;
      double ddz = (dpz - dnz) / grid->cell_size;
      // prevent negative concentrations
      grid->conc[isub][i] = MAX(THRESHOLD_CONC, prev[i] + dt * (ddx + ddy + ddz + grid->reac[isub][i]));
    }
  }
}

/* ----------------------------------------------------------------------
 Average substrate distribution before solving diffusion in closed system.
 ------------------------------------------------------------------------- */
void FixDiffusionReaction::closed_system_initial()
{
  if (!closed_system) return;
  double ave_conc = 0;
  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GRID_MASK) {
      ave_conc += grid->conc[isub][i];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &ave_conc, 1, MPI_DOUBLE, MPI_SUM, world);
  ave_conc /= (grid->box[0] * grid->box[1] * grid->box[2]);
  for (int i = 0; i < grid->ncells; i++) {
    grid->conc[isub][i] = ave_conc;
  }
}

/* ----------------------------------------------------------------------
 Scaleup substrate concentrations in closed system based on residual value and
 biological timestep
 ------------------------------------------------------------------------- */
void FixDiffusionReaction::closed_system_scaleup(double biodt)
{
  if (!closed_system) return;
  for (int i = 0; i < grid->ncells; i++) {
    double res = grid->conc[isub][i] - prev[i];
    grid->conc[isub][i] += res / dt * biodt;
    grid->conc[isub][i] = MAX(0, grid->conc[isub][i]);
  }
}
