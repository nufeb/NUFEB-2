/*
 * fix_kinetics/diffusionS.cpp
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include "fix_bio_kinetics_diffusion.h"

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "atom.h"
#include "domain.h"
#include "error.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "atom_vec_bio.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "memory.h"
#include "update.h"
#include "variable.h"
#include "group.h"
#include "comm.h"

#define BUFMIN 1000

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion::FixKineticsDiffusion(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {

  if (narg < 8)
    error->all(FLERR, "Not enough arguments in fix diffusion command");

  shearflag = dragflag = 0;
  bulkflag = 0;
  srate = 0;

  var = new char*[1];
  ivar = new int[1];

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[3 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3 + i][2]);
  }

  //set boundary condition flag:
  //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  if (strcmp(arg[4], "pp") == 0)
    xbcflag = 0;
  else if (strcmp(arg[4], "dd") == 0)
    xbcflag = 1;
  else if (strcmp(arg[4], "nd") == 0)
    xbcflag = 2;
  else if (strcmp(arg[4], "nn") == 0)
    xbcflag = 3;
  else if (strcmp(arg[4], "dn") == 0)
    xbcflag = 4;
  else
    error->all(FLERR, "Illegal x-axis boundary condition command");

  if (strcmp(arg[5], "pp") == 0)
    ybcflag = 0;
  else if (strcmp(arg[5], "dd") == 0)
    ybcflag = 1;
  else if (strcmp(arg[5], "nd") == 0)
    ybcflag = 2;
  else if (strcmp(arg[5], "nn") == 0)
    ybcflag = 3;
  else if (strcmp(arg[5], "dn") == 0)
    ybcflag = 4;
  else
    error->all(FLERR, "Illegal y-axis boundary condition command");

  if (strcmp(arg[6], "pp") == 0)
    zbcflag = 0;
  else if (strcmp(arg[6], "dd") == 0)
    zbcflag = 1;
  else if (strcmp(arg[6], "nd") == 0)
    zbcflag = 2;
  else if (strcmp(arg[6], "nn") == 0)
    zbcflag = 3;
  else if (strcmp(arg[6], "dn") == 0)
    zbcflag = 4;
  else
    error->all(FLERR, "Illegal z-axis boundary condition command");

  if (strcmp(arg[7], "kg") == 0)
    unit = 1;
  else if (strcmp(arg[7], "mol") == 0)
    unit = 0;
  else
    error->all(FLERR, "Illegal unit in fix kinetics/diffusion command: specify 'kg' or 'mol'");

  int iarg = 8;
  while (iarg < narg){
    if (strcmp(arg[iarg],"srate") == 0) {
      srate = force->numeric(FLERR, arg[iarg+1]);
      if (srate < 0)
        error->all(FLERR, "Illegal fix kinetics/diffusion command: srate");
      iarg += 2;
    } else if (strcmp(arg[iarg],"bulk") == 0) {
      bulkflag = 1;
      q = force->numeric(FLERR, arg[iarg+1]);
      if (q < 0)
        error->all(FLERR, "Flow rate (frate) cannot be negative");
      rvol = force->numeric(FLERR, arg[iarg+2]);
      if (rvol <= 0)
        error->all(FLERR, "Reactor volume (rvol) cannot be equal or less than 0");
      af = force->numeric(FLERR, arg[iarg+3]);
      if (af < 0)
        lmp->error->all(FLERR, "Biofilm surface area (Af) cannot be negative");
      iarg += 4;
    }
    else
      error->all(FLERR, "Illegal fix kinetics/diffusion command");
  }
}

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion::~FixKineticsDiffusion() {
  for (int i = 0; i < 1; i++) {
    delete[] var[i];
  }

  delete[] var;
  delete[] ivar;

  memory->destroy(diffD);
  memory->destroy(xGrid);
  memory->destroy(nuGrid);
  memory->destroy(nuPrev);
  memory->destroy(ghost);

  delete [] requests;
}

/* ---------------------------------------------------------------------- */

int FixKineticsDiffusion::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsDiffusion::init() {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  if (!atom->radius_flag)
    error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix diffusion does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix diffusion is invalid style");
  }

  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "shear") == 0) {
      shearflag = 1;
    } else if (strcmp(modify->fix[j]->style, "fdrag") == 0) {
      dragflag = 1;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR, "The fix kinetics command is required for running iBM simulation");

  bio = kinetics->bio;
  tol = input->variable->compute_equal(ivar[0]);

  //set diffusion grid size
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;
  bzhi = kinetics->bnz * stepz;

  if (!isEuqal(stepx, stepy, stepz))
    error->all(FLERR, "Grid is not cubic");

  nX = kinetics->subn[0] + 2;
  nY = kinetics->subn[1] + 2;
  nZ = kinetics->subn[2] + 2;

  nXYZ = nX * nY * nZ;

  int nnus = bio->nnus;
  diffD = memory->create(diffD, nnus + 1, "diffusion:diffD");
  //inlet concentration, diffusion constant
  //and maximum boundary condition conc value
  for (int i = 1; i <= nnus; i++) {
    diffD[i] = bio->diffCoeff[i];
  }

  xGrid = memory->create(xGrid, nXYZ, 3, "diffusion:xGrid");
  nuGrid = memory->create(nuGrid, nnus + 1, nXYZ, "diffusion:nuGrid");
  nuPrev = memory->create(nuPrev, nnus + 1, nXYZ, "diffusion:nuPrev");
  ghost = memory->create(ghost, nXYZ, "diffusion:ghost");

  init_grid();

  // create request vector
  requests = new MPI_Request[MAX(2 * comm->nprocs, nnus + 1)];

  setup_exchange(kinetics->grid, kinetics->subgrid.get_box(), {xbcflag == 0, ybcflag == 0, zbcflag == 0});
}

/* ----------------------------------------------------------------------
 solve diffusion and reaction
 ------------------------------------------------------------------------- */

int *FixKineticsDiffusion::diffusion(int *nuConv, int iter, double diffT) {
  if (iter == 1 && kinetics->blayer >= 0)
    update_grids();

  int nnus = bio->nnus;
  this->diffT = diffT;
  double **nuR = kinetics->nuR;
  double **nuS = kinetics->nuS;
  double *nuBS = kinetics->nuBS;
  double **iniS = bio->iniS;

  DecompGrid<FixKineticsDiffusion>::exchange();

  double *maxS = new double[nnus + 1];
  for (int i = 1; i <= nnus; i++) {
    if (bio->nuType[i] == 0 && !nuConv[i]) {
      if (unit == 0) {
        xbcm = iniS[i][1] * 1000;
        xbcp = iniS[i][2] * 1000;
        ybcm = iniS[i][3] * 1000;
        ybcp = iniS[i][4] * 1000;
        zbcm = iniS[i][5] * 1000;
        zbcp = iniS[i][6] * 1000;
      } else {
        xbcm = iniS[i][1];
        xbcp = iniS[i][2];
        ybcm = iniS[i][3];
        ybcp = iniS[i][4];
        zbcm = iniS[i][5];
        zbcp = iniS[i][6];
      }
      // copy current concentrations
      for (int grid = 0; grid < nXYZ; grid++) {
        nuPrev[i][grid] = nuGrid[i][grid];
      }

      maxS[i] = 0;
      // solve diffusion and reaction
      for (int grid = 0; grid < nXYZ; grid++) {
        // transform nXYZ index to nuR index
        if (!ghost[grid]) {
          int ind = get_index(grid);
          double r = (unit == 1) ? (r = nuR[i][ind]) :  (r = nuR[i][ind] * 1000);

          compute_flux(diffD[i], nuGrid[i][grid], nuPrev[i], r, grid, ind);

          nuR[i][ind] = 0;

          if (nuGrid[i][grid] > 0) {
            (unit == 1) ? (nuS[i][ind] = nuGrid[i][grid]) : (nuS[i][ind] = nuGrid[i][grid] / 1000);
          } else {
            nuGrid[i][grid] = 1e-20;
            nuS[i][ind] = 1e-20;
          }

	  if (maxS[i] < nuGrid[i][grid])
	    maxS[i] = nuGrid[i][grid];
        } else
          compute_bc(nuGrid[i][grid], nuPrev[i], grid, nuBS[i]);
      }
#if MPI_VERSION >= 3
      MPI_Iallreduce(MPI_IN_PLACE, &maxS[i], 1, MPI_DOUBLE, MPI_MAX, world, &requests[i]);
#else
      MPI_Allreduce(MPI_IN_PLACE, &maxS[i], 1, MPI_DOUBLE, MPI_MAX, world);
#endif
    }
  }

  int nrequests = 0;
  for (int i = 1; i <= nnus; i++) {
    if (bio->nuType[i] == 0 && !nuConv[i]) { // checking if is liquid
#if MPI_VERSION >= 3
      MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
#endif
      // check convergence criteria
      nuConv[i] = true;
      for (int grid = 0; grid < nXYZ; grid++) {
        if (!ghost[grid]) {
          double ratio = 1000;
          double div = (maxS[i] == 0) ? 1 : maxS[i];
          double rate = nuGrid[i][grid] / div;
          double prevRate = nuPrev[i][grid] / div;
          ratio = fabs(rate - prevRate);

          nuConv[i] &= ratio < tol;
          if (!nuConv[i])
            break;
        }
      }
#if MPI_VERSION >= 3
      MPI_Iallreduce(MPI_IN_PLACE, &nuConv[i], 1, MPI_INT, MPI_BAND, world, &requests[nrequests++]);
#else
      MPI_Allreduce(MPI_IN_PLACE, &nuConv[i], 1, MPI_INT, MPI_BAND, world);
#endif
    }
  }
#if MPI_VERSION >= 3
  MPI_Waitall(nrequests, requests, MPI_STATUS_IGNORE);
#endif

  delete [] maxS;

  return nuConv;
}

/* ----------------------------------------------------------------------
 Update attributes in non-boundary grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::update_grids() {

  bzhi = kinetics->bnz * stepz;
  if (kinetics->bgrids == 0)
    nXYZ = 0;
  else
    nXYZ = nX * nY * (MIN(kinetics->subn[2], MAX(0, kinetics->bnz - kinetics->subnlo[2])) + 2);

  for (int grid = 0; grid < nXYZ; grid++) {
    if (xGrid[grid][0] < kinetics->sublo[0]|| xGrid[grid][1] < kinetics->sublo[1] || xGrid[grid][2] < kinetics->sublo[2]
    || xGrid[grid][0] > kinetics->subhi[0] || xGrid[grid][1] > kinetics->subhi[1] || xGrid[grid][2] > MIN(bzhi, kinetics->subhi[2]))
      ghost[grid] = true;
    else
      ghost[grid] = false;
  }
}

/* ----------------------------------------------------------------------
 Initialize grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::init_grid() {
  double *nuBS = kinetics->nuBS;
  double **iniS = bio->iniS;

  double i, j, k;
  int cell = 0;
  for (k = kinetics->sublo[2] - (stepz/2); k < kinetics->subhi[2] + stepz; k += stepz) {
    for (j = kinetics->sublo[1] - (stepy/2); j < kinetics->subhi[1] + stepy; j += stepy) {
      for (i = kinetics->sublo[0] - (stepx/2); i < kinetics->subhi[0] + stepx; i += stepx) {
        xGrid[cell][0] = i;
        xGrid[cell][1] = j;
        xGrid[cell][2] = k;
        //Initialise concentration values for ghost and std grids
        for (int nu = 1; nu <= bio->nnus; nu++) {
          if (i < kinetics->sublo[0]) {
            ghost[cell] = true;
            nuGrid[nu][cell] = iniS[nu][1];
          } else if (i > kinetics->subhi[0]) {
            ghost[cell] = true;
            nuGrid[nu][cell] = iniS[nu][2];
          } else if (j < kinetics->sublo[1]) {
            ghost[cell] = true;
            nuGrid[nu][cell] = iniS[nu][3];
          } else if (j > kinetics->subhi[1]) {
            ghost[cell] = true;
            nuGrid[nu][cell] = iniS[nu][4];
          } else if (k < kinetics->sublo[2]) {
            ghost[cell] = true;
            nuGrid[nu][cell] = iniS[nu][5];
          } else if (k > MIN(bzhi, kinetics->subhi[2])) {
            ghost[cell] = true;
            nuGrid[nu][cell] = iniS[nu][6];
          } else {
	    ghost[cell] = false;
	    nuGrid[nu][cell] = iniS[nu][0];
          }

	  if (unit == 0)
	    nuGrid[nu][cell] = nuGrid[nu][cell] * 1000;
	  if (cell == 0 && unit == 1)
	    nuBS[nu] = iniS[nu][6];
	  else if (cell == 0 && unit == 0)
	    nuBS[nu] = iniS[nu][6] * 1000;
        }
        cell++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 Mass balance in bulk liquid
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bulk() {
  double vol = stepx * stepy * stepz;
  double **iniS = bio->iniS;
  double **nuR = kinetics->nuR;
  double *nuBS = kinetics->nuBS;

  for (int nu = 1; nu <= bio->nnus; nu++) {
    // the concentration of o2 in the bulk liquid is kept constant by aeration
    if (strcmp(bio->nuName[nu], "o2") != 0) {
      double sumR = 0;
      double global_sumR = 0;
      double iniBC = (unit == 1) ? iniBC = iniS[nu][6] : iniBC = iniS[nu][6] * 1000;
      // sum up consumption
      for (int i = 0; i < kinetics->bgrids; i++) {
        (unit == 1) ? (sumR += nuR[nu][i]) : (sumR += nuR[nu][i] * 1000);
      }

      MPI_Allreduce(&sumR, &global_sumR, 1, MPI_DOUBLE, MPI_SUM, world);

      double dt = update->dt * kinetics->nevery;
      // solve for the mass balance in bulk liquid

      nuBS[nu] = nuBS[nu] + ((q / rvol) * (iniBC - nuBS[nu]) + ((af * global_sumR * vol) / (rvol * yhi * xhi))) * dt;

      if (nuBS[nu] < 0) nuBS[nu] = 1e-20;
    }
  }
}

/* ----------------------------------------------------------------------
 get index of non-ghost mesh grid
 ------------------------------------------------------------------------- */

int FixKineticsDiffusion::get_index(int grid) {
  int ix = (xGrid[grid][0] - kinetics->sublo[0]) / stepx;
  int iy = (xGrid[grid][1] - kinetics->sublo[1]) / stepy;
  int iz = (xGrid[grid][2] - kinetics->sublo[2]) / stepz;

  int ind = iz * kinetics->subn[0] * kinetics->subn[1] + iy * kinetics->subn[0] + ix;

  return ind;
}

/* ----------------------------------------------------------------------
 update concentration for ghost grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bc(double &nuCell, double *nuPrev, int grid, double bulk) {
  //for nx = ny = nz = 1 grids
  //18 19 20        21 22 23       24 25 26
  //9  10 11        12 13 14       15 16 17
  //0   1  2         3  4  5        6  7  8
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - nX;  // y direction
  int fwd = grid + nX;  // y direction
  int down = grid - nX * nY; // z direction
  int up = grid + nX * nY;  // z dirction

  // assign values to the ghost-grids according to the boundary conditions.
  // If ghostcells are Neu then take the values equal from the adjacent cells.
  // if ghostcells are dirich then take the values equal to negative of the adjacent cells.
  // if ghostcells are mixed then zlo ghost cells are nuemann, zhi ghost cells are dirichlet, other four surfaces are periodic BC.
  // low-z surface
  if (xGrid[grid][2] < zlo && !ghost[up]) {
    //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
    if (zbcflag == 0 && kinetics->nz == kinetics->subn[2]) {
      int zhiGrid = grid + nX * nY * nz;
      nuCell = nuPrev[zhiGrid];
    } else if (zbcflag == 1) {
      nuCell = 2 * zbcm - nuPrev[up];
    } else if (zbcflag == 2) {
      nuCell = nuPrev[up];
    } else if (zbcflag == 3) {
      nuCell = nuPrev[up];
    } else if (zbcflag == 4) {
      nuCell = 2 * zbcm - nuPrev[up];
    }
  }
  // high-z surface
  else if (xGrid[grid][2] > bzhi && !ghost[down]) {
    if (zbcflag == 0 && kinetics->nz == kinetics->subn[2]) {
      int zloGrid = grid - nX * nY * nz;
      nuCell = nuPrev[zloGrid];
    } else if (zbcflag == 1) {
      nuCell = 2 * bulk - nuPrev[down];
    } else if (zbcflag == 2) {
      nuCell = 2 * bulk - nuPrev[down];
    } else if (zbcflag == 3) {
      nuCell = nuPrev[down];
    } else if (zbcflag == 4) {
      nuCell = nuPrev[down];
    }
  }
  // low-y surface
  else if (xGrid[grid][1] < ylo && !ghost[fwd]) {
    if (ybcflag == 0 && kinetics->ny == kinetics->subn[1]) {
      int yhiGrid = grid + nX * ny;
      nuCell = nuPrev[yhiGrid];
    } else if (ybcflag == 1) {
      nuCell = 2 * ybcm - nuPrev[fwd];
    } else if (ybcflag == 2) {
      nuCell = nuPrev[fwd];
    } else if (ybcflag == 3) {
      nuCell = nuPrev[fwd];
    } else if (ybcflag == 4) {
      nuCell = 2 * ybcm - nuPrev[fwd];
    }
  }
  // high-y surface
  else if (xGrid[grid][1] > yhi && !ghost[bwd]) {
    if (ybcflag == 0 && kinetics->ny == kinetics->subn[1]) {
      int yloGrid = grid - nX * ny;
      nuCell = nuPrev[yloGrid];
    } else if (ybcflag == 1) {
      nuCell = 2 * ybcp - nuPrev[bwd];
    } else if (ybcflag == 2) {
      nuCell = 2 * ybcp - nuPrev[bwd];
    } else if (ybcflag == 3) {
      nuCell = nuPrev[bwd];
    } else if (ybcflag == 4) {
      nuCell = nuPrev[bwd];
    }
  }
  // low-x surface
  else if (xGrid[grid][0] < xlo && !ghost[rhs]) {
    if (xbcflag == 0 && kinetics->nx == kinetics->subn[0]) {
      int xhiGrid = grid + nx;
      nuCell = nuPrev[xhiGrid];
    } else if (xbcflag == 1) {
      nuCell = 2 * xbcm - nuPrev[rhs];
    } else if (xbcflag == 2) {
      nuCell = nuPrev[rhs];
    } else if (xbcflag == 3) {
      nuCell = nuPrev[rhs];
    } else if (xbcflag == 4) {
      nuCell = 2 * xbcm - nuPrev[rhs];
    }
  }
  // high-x surface
  else if (xGrid[grid][0] > xhi && !ghost[lhs]) {
    if (xbcflag == 0 && kinetics->nx == kinetics->subn[0]) {
      int xloGrid = grid - nx;
      nuCell = nuPrev[xloGrid];
    } else if (xbcflag == 1) {
      nuCell = 2 * xbcp - nuPrev[lhs];
    } else if (xbcflag == 2) {
      nuCell = 2 * xbcp - nuPrev[lhs];
    } else if (xbcflag == 3) {
      nuCell = nuPrev[lhs];
    } else if (xbcflag == 4) {
      nuCell = nuPrev[lhs];
    }
  }
}

/* ----------------------------------------------------------------------
 update concentration for non-ghost grids
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_flux(double cellDNu, double &nuCell, double *nuPrev, double rateNu, int grid, int ind) {
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - nX;  // y direction
  int fwd = grid + nX;  // y direction
  int down = grid - nX * nY; // z direction
  int up = grid + nX * nY;  // z direction

  double jRight = cellDNu * (nuPrev[rhs] - nuPrev[grid]) / stepx;
  double jLeft = cellDNu * (nuPrev[grid] - nuPrev[lhs]) / stepx;
  double jX = (jRight - jLeft) / stepx;

  double jForward = cellDNu * (nuPrev[fwd] - nuPrev[grid]) / stepy;
  double jBackward = cellDNu * (nuPrev[grid] - nuPrev[bwd]) / stepy;
  double jY = (jForward - jBackward) / stepy;

  double jUp = cellDNu * (nuPrev[up] - nuPrev[grid]) / stepz;
  double jDown = cellDNu * (nuPrev[grid] - nuPrev[down]) / stepz;
  double jZ = (jUp - jDown) / stepz;

  double res = 0;
  double shear = 0;

  // Adding fluxes in all the directions and the uptake rate (RHS side of the equation)

  if (dragflag == 1) {
    double uX = kinetics->fV[0][ind] * (nuPrev[rhs] - nuPrev[lhs]) / stepx;
    double uY = kinetics->fV[1][ind] * (nuPrev[fwd] - nuPrev[bwd]) / stepy;
    double uZ = kinetics->fV[2][ind] * (nuPrev[up] - nuPrev[down]) / stepz;

    res = (jX + jY + jZ + rateNu - uX - uY - uZ) * diffT;
   // printf("grid = %i, ind = %i, %e %e %e \n", grid, ind, kinetics->fV[0][ind], kinetics->fV[1][ind], kinetics->fV[2][ind]);
  } else if (shearflag == 1) {
    int hgrid = grid;
    int deep = 0;

    while (!ghost[hgrid]) {
      hgrid = hgrid - nX * nY;
      deep++;
    }

    shear = srate * (deep * stepz - stepz / 2) * (nuPrev[rhs] - nuPrev[lhs]) / (2 * stepz);

    res = (jX + jY + jZ + rateNu - shear) * diffT;
  } else {
    res = (jX + jY + jZ + rateNu) * diffT;
  }

  //Updating the value: Ratesub*diffT + nuCell[cell](previous)
  nuCell = nuPrev[grid] + res;
}

/* ----------------------------------------------------------------------
 Compare double values for equality
 ------------------------------------------------------------------------- */

bool FixKineticsDiffusion::isEuqal(double a, double b, double c) {
  double epsilon = 1e-10;
  if ((fabs(a - b) > epsilon) || (fabs(a - c) > epsilon) || (fabs(b - c) > epsilon))
    return false;

  return true;
}

/* ----------------------------------------------------------------------
 Manually update reaction if none of the surface is using dirichlet bc
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::update_nuS() {
  if ((xbcflag == 0 || xbcflag == 3) && (ybcflag == 0 || ybcflag == 3) && (zbcflag == 0 || zbcflag == 3)) {
    double **nuS = kinetics->nuS;
    double **nuR = kinetics->nuR;

    for (int nu = 1; nu < bio->nnus + 1; nu++) {
      for (int grid = 0; grid < nXYZ; grid++) {
        // transform nXYZ index to nuR index
        if (!ghost[grid]) {
          int ind = get_index(grid);
          int dt = update->dt * kinetics->nevery;
          double r = (unit == 1) ? nuR[nu][ind] * dt : nuR[nu][ind] *dt * 1000;

          nuGrid[nu][grid] += r;

          if (nuGrid[nu][grid] <= 0)
            nuGrid[nu][grid] = 1e-20;

	  (unit == 1) ? nuS[nu][ind] = nuGrid[nu][grid] : nuS[nu][ind] = nuGrid[nu][grid] / 1000;
        }
      }
    }
  }
}

int FixKineticsDiffusion::get_elem_per_cell() const {
  return bio->nnus;
}

void FixKineticsDiffusion::resize(const Subgrid<double, 3> &subgrid) {
  int nnus = bio->nnus;
  nXYZ = subgrid.cell_count();
  xGrid = memory->grow(xGrid, nXYZ, 3, "diffusion:xGrid");
  nuGrid = memory->grow(nuGrid, nnus + 1, nXYZ, "diffusion:nuGrid");
  nuPrev = memory->grow(nuPrev, nnus + 1, nXYZ, "diffusion:nuPrev");
  ghost = memory->grow(ghost, nXYZ, "diffusion:ghost");
}

void FixKineticsDiffusion::migrate(const Grid<double, 3> &grid, const Box<int, 3> &from, const Box<int, 3> &to) {
  DecompGrid<FixKineticsDiffusion>::migrate(grid, from, to, extend(from), extend(to));
  Subgrid<double, 3> subgrid(grid, to);
  setup_exchange(grid, to, {xbcflag == 0, ybcflag == 0, zbcflag == 0});
  Subgrid<double, 3> extended(kinetics->grid, extend(to));
  auto cell_centers = extended.get_cell_centers();
  for (int i = 0; i < extended.cell_count(); i++) {
    xGrid[i][0] = cell_centers[i][0];
    xGrid[i][1] = cell_centers[i][1];
    xGrid[i][2] = cell_centers[i][2];
  }
  nx = subgrid.get_dimensions()[0];
  ny = subgrid.get_dimensions()[1];
  nz = subgrid.get_dimensions()[2];
  nX = extended.get_dimensions()[0];
  nY = extended.get_dimensions()[1];
  nZ = extended.get_dimensions()[2];
  for (int grid = 0; grid < nXYZ; grid++) {
    if (xGrid[grid][0] < kinetics->sublo[0] || xGrid[grid][1] < kinetics->sublo[1] || xGrid[grid][2] < kinetics->sublo[2] ||
	xGrid[grid][0] > kinetics->subhi[0] || xGrid[grid][1] > kinetics->subhi[1] || xGrid[grid][2] > kinetics->subhi[2])
      ghost[grid] = true;
    else
      ghost[grid] = false;
  }
  double *nuBS = kinetics->nuBS;
  for (int i = 1; i <= bio->nnus; i++) {
    for (int grid = 0; grid < nXYZ; grid++) {
      if (ghost[grid]) {
	compute_bc(nuGrid[i][grid], nuGrid[i], grid, nuBS[i]);
      }
    }
  }
}
