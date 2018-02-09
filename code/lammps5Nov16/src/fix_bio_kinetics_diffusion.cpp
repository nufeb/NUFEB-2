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

#include "thr_omp.h"

#define BUFMIN 1000

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion::FixKineticsDiffusion(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {

  if (narg != 12)
    error->all(FLERR, "Not enough arguments in fix diffusion command");

  var = new char*[5];
  ivar = new int[5];

  for (int i = 0; i < 5; i++) {
    int n = strlen(&arg[3 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[3 + i][2]);
  }

  //set boundary condition flag:
  //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  if (strcmp(arg[8], "pp") == 0)
    xbcflag = 0;
  else if (strcmp(arg[8], "dd") == 0)
    xbcflag = 1;
  else if (strcmp(arg[8], "nd") == 0)
    xbcflag = 2;
  else if (strcmp(arg[8], "nn") == 0)
    xbcflag = 3;
  else if (strcmp(arg[8], "dn") == 0)
    xbcflag = 4;
  else
    error->all(FLERR, "Illegal x-axis boundary condition command");

  if (strcmp(arg[9], "pp") == 0)
    ybcflag = 0;
  else if (strcmp(arg[9], "dd") == 0)
    ybcflag = 1;
  else if (strcmp(arg[9], "nd") == 0)
    ybcflag = 2;
  else if (strcmp(arg[9], "nn") == 0)
    ybcflag = 3;
  else if (strcmp(arg[9], "dn") == 0)
    ybcflag = 4;
  else
    error->all(FLERR, "Illegal y-axis boundary condition command");

  if (strcmp(arg[10], "pp") == 0)
    zbcflag = 0;
  else if (strcmp(arg[10], "dd") == 0)
    zbcflag = 1;
  else if (strcmp(arg[10], "nd") == 0)
    zbcflag = 2;
  else if (strcmp(arg[10], "nn") == 0)
    zbcflag = 3;
  else if (strcmp(arg[10], "dn") == 0)
    zbcflag = 4;
  else if (strcmp(arg[10], "db") == 0)
    zbcflag = 5;
  else
    error->all(FLERR, "Illegal z-axis boundary condition command");

  if (strcmp(arg[11], "kg") == 0)
    unit = 1;
  else if (strcmp(arg[11], "mol") == 0)
    unit = 0;
  else
    error->all(FLERR, "Illegal unit in fix kinetics/diffusionS command: specify 'kg' or 'mol'");

  shearflag = dragflag = 0;

//  rstep = atoi(arg[11]);
//  if (rstep > 1) rflag = 1;
}

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion::~FixKineticsDiffusion() {
  int i;
  for (i = 0; i < 5; i++) {
    delete[] var[i];
  }

  delete[] var;
  delete[] ivar;

  memory->destroy(diffD);
  memory->destroy(xGrid);
  memory->destroy(nuGrid);
  memory->destroy(nuPrev);
  memory->destroy(ghost);
  memory->destroy(nuBS);

  memory->destroy(recvbuff);
  memory->destroy(sendbuff);
  memory->destroy(convergences);

  delete [] requests;
  delete [] status;
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

  for (int n = 0; n < 5; n++) {
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
  shearRate = input->variable->compute_equal(ivar[0]);
  tol = input->variable->compute_equal(ivar[1]);
  q = input->variable->compute_equal(ivar[2]);
  rvol = input->variable->compute_equal(ivar[3]);
  if (rvol <= 0) lmp->error->all(FLERR, "Reactor volume cannot be equal or less than 0");
  af = input->variable->compute_equal(ivar[4]);

  //set diffusion grid size
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  nnus = bio->nnus;
  iniS = bio->iniS;
  diffCoeff = bio->diffCoeff;

  diffD = memory->create(diffD, nnus + 1, "diffusion:diffD");

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

  //inlet concentration, diffusion constant
  //and maximum boundary condition conc value
  for (int i = 1; i <= nnus; i++) {
    diffD[i] = diffCoeff[i];
  }

  xGrid =  memory->create(xGrid, nXYZ, 3, "diffusion:xGrid");
  nuGrid = memory->create(nuGrid, nnus + 1, nXYZ, "diffusion:nuGrid");
  nuPrev = memory->create(nuPrev, nXYZ, "diffusion:nuPrev");
  ghost = memory->create(ghost, nXYZ, "diffusion:ghost");
  nuBS = memory->create(nuBS, nnus + 1, "diffusion:nuBS");

  // TODO: optimize for first-touch policy
  //initialise grids
  double i, j, k;
  int grid = 0;
  for (k = kinetics->sublo[2] - (stepz/2); k < kinetics->subhi[2] + stepz; k += stepz) {
    for (j = kinetics->sublo[1] - (stepy/2); j < kinetics->subhi[1] + stepy; j += stepy) {
      for (i = kinetics->sublo[0] - (stepx/2); i < kinetics->subhi[0] + stepx; i += stepx) {
        xGrid[grid][0] = i;
        xGrid[grid][1] = j;
        xGrid[grid][2] = k;
        //Initialise concentration values for ghost and std grids
        for (int nu = 1; nu <= nnus; nu++) {
          if (i < kinetics->sublo[0]) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][1];
          } else if (i > kinetics->subhi[0]) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][2];
          } else if (j < kinetics->sublo[1]) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][3];
          } else if (j > kinetics->subhi[1]) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][4];
          } else if (k < kinetics->sublo[2]) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][5];
          } else if (k > bzhi) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][6];
          } else {
            ghost[grid] = false;
            nuGrid[nu][grid] = iniS[nu][0];
          }

          if (unit == 0)
            nuGrid[grid][nu] = nuGrid[grid][nu] * 1000;
          if (grid == 0 && unit == 1)
            nuBS[nu] = iniS[nu][6];
          else if (grid == 0 && unit == 0)
            nuBS[nu] = iniS[nu][6] * 1000;
        }
        grid++;
      }
    }
  }

  // create send and recv buffers
  recv_buff_size = BUFMIN;
  recvbuff = memory->create(recvbuff, recv_buff_size, "diffusion::recvbuff");
  send_buff_size = BUFMIN;
  sendbuff = memory->create(sendbuff, send_buff_size, "diffusion::sendbuff");
  convergences = memory->create(convergences, comm->nprocs, "diffusion::convergences");

  // create request vector
  requests = new MPI_Request[MAX(comm->nprocs, nnus + 1)];
  status = new MPI_Status[MAX(comm->nprocs, nnus + 1)];
}

/* ----------------------------------------------------------------------
 solve diffusion and reaction
 ------------------------------------------------------------------------- */

bool* FixKineticsDiffusion::diffusion(bool *nuConv, int iter, double diffT) {
  if (iter == 1 && kinetics->bl > 0)
    update_grids();

  this->diffT = diffT;
  nuS = kinetics->nuS;
  nuR = kinetics->nuR;

  int nrecvcells = kinetics->recvend[comm->nprocs - 1];
  int nsendcells = kinetics->sendend[comm->nprocs - 1];
  if (recv_buff_size < nrecvcells) {
    recv_buff_size += ((nrecvcells * nnus) / BUFMIN + 1) * BUFMIN;
    memory->grow(recvbuff, recv_buff_size, "diffusion::recvbuff");
  }
  if (send_buff_size < nsendcells) {
    send_buff_size += ((nsendcells * nnus) / BUFMIN + 1) * BUFMIN;
    memory->grow(sendbuff, send_buff_size, "diffusion::recvbuff");
  }

  // copy nutrient grid data to send buffer
  for (int c = 0; c < nsendcells; c++) {
    int cell = kinetics->sendcells[c];
    for (int n = 0; n < nnus; n++) {
      sendbuff[c * nnus + n] = nuGrid[n][cell]; 
    }  
  }
  // send and recv grid data
  int nrequests = 0;
  for (int p = 0; p < comm->nprocs; p++) {
    if (p == comm->me)
      continue;
    int recvn = (kinetics->recvend[p] - kinetics->recvbegin[p]) * nnus;
    if (recvn > 0)
      MPI_Irecv(&recvbuff[kinetics->recvbegin[p] * nnus], recvn, MPI_DOUBLE, p, 0, world, &requests[nrequests++]);
    int sendn = (kinetics->sendend[p] - kinetics->sendbegin[p]) * nnus;
    if (sendn > 0)
      MPI_Isend(&sendbuff[kinetics->sendbegin[p] * nnus], sendn, MPI_DOUBLE, p, 0, world, &requests[nrequests++]);
  }
  // wait for all MPI requests
  MPI_Waitall(nrequests, requests, status);
  // copy received data to nuGrid
  for (int c = 0; c < nrecvcells; c++) {
    int cell = kinetics->recvcells[c];
    for (int n = 0; n < nnus; n++) {
      nuGrid[n][cell] = recvbuff[c * nnus + n];
    }
  }

  for (int i = 1; i <= nnus; i++) {
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
    
    if(iter == 1 && strcmp(bio->nuName[i], "o2") != 0 && q >= 0 && af >= 0) compute_bulk(i);
  }

  double *maxS = new double[nnus + 1];
  for (int i = 1; i <= nnus; i++) {
    if (bio->nuType[i] == 0 && !nuConv[i]) { // checking if is liquid
      // copy current concentrations
      for (int grid = 0; grid < nXYZ; grid++) {
	nuPrev[grid] = nuGrid[i][grid];
      }
	
      maxS[i] = 0;
      // solve diffusion and reaction
      for (int grid = 0; grid < nXYZ; grid++) {
	// transform nXYZ index to nuR index
	if (!ghost[grid]) {
	  int ix = floor((xGrid[grid][0] - kinetics->sublo[0])/stepx);
	  int iy = floor((xGrid[grid][1] - kinetics->sublo[1])/stepy);
	  int iz = floor((xGrid[grid][2] - kinetics->sublo[2])/stepz);
	  int ind = iz * kinetics->subn[0] * kinetics->subn[1] + iy * kinetics->subn[0] + ix;

          double r;
          if (unit == 0) r = nuR[i][ind] * 1000;
          else r = nuR[i][ind];

	  compute_flux(diffD[i], nuGrid[i][grid], nuPrev, r, grid);

	  nuR[i][ind] = 0;

	  if (nuGrid[i][grid] > 0) {
	    if (unit == 0) nuS[i][ind] = nuGrid[i][grid] / 1000;
	    else nuS[i][ind] = nuGrid[i][grid];
	  }
	  else {
	    nuGrid[i][grid] = 1e-20;
	    nuS[i][ind] = 1e-20;
	  }
	} else compute_bc(nuGrid[i][grid], nuPrev, grid, nuBS[i]);

	if (maxS[i] < nuGrid[i][grid]) maxS[i] = nuGrid[i][grid];
      }
      MPI_Iallreduce(MPI_IN_PLACE, &maxS[i], 1, MPI_DOUBLE, MPI_MAX, world, &requests[i]);
    }
  }

  nrequests = 0;
  for (int i = 1; i <= nnus; i++) {
    if (bio->nuType[i] == 0 && !nuConv[i]) { // checking if is liquid
      MPI_Wait(&requests[i], &status[i]);
      // check convergence criteria
      nuConv[i] = true;
      for (int grid = 0; grid < nXYZ; grid++) {
	if(!ghost[grid]) {
	  double ratio = 1000;
	  double div = (maxS[i] == 0) ? 1 : maxS[i]; 
	  double rate = nuGrid[i][grid] / div;
	  double prevRate = nuPrev[grid] / div;
	  ratio = fabs(rate - prevRate);

	  nuConv[i] &= (ratio < tol);
	  if (!nuConv[i]) break;
	}
      }
      MPI_Iallreduce(MPI_IN_PLACE, &nuConv[i], 1, MPI_INT, MPI_BAND, world, &requests[nrequests++]);
    }
  }
  MPI_Waitall(nrequests, requests, status);

  delete [] maxS;

  return nuConv;
}

void FixKineticsDiffusion::update_grids() {
  //update grids
  bzhi = kinetics->bnz * stepz;
  nXYZ = nX * nY * (kinetics->bnz + 2);

  for (int grid = 0; grid < nXYZ; grid++) {
    if (xGrid[grid][0] < 0 || xGrid[grid][1] < 0 || xGrid[grid][2] < 0
        || xGrid[grid][0] > kinetics->subhi[0] || xGrid[grid][1] > kinetics->subhi[1] || xGrid[grid][2] > bzhi)
      ghost[grid] = true;
    else
      ghost[grid] = false;
  }
}

/* ----------------------------------------------------------------------
 Mass balances of nutrients in the bulk liquid
 ------------------------------------------------------------------------- */

void FixKineticsDiffusion::compute_bulk(int nu) {
  double sumR = 0;
  double vol = stepx * stepy * stepz;

  for (int i = 0; i < nx * ny * kinetics->nz; i++) {
    if (unit == 0) sumR += nuR[nu][i] * vol * 1000;
    else sumR += nuR[nu][i] * vol;
  }

  nuBS[nu] = nuBS[nu] + ((q / rvol) * (zbcp - nuBS[nu]) + (af / (rvol * yhi * xhi)) * sumR) * update->dt * nevery;
  printf("bulk liquid = %e \n", nuBS[nu]);
}

/* ----------------------------------------------------------------------
  get index of non-ghost mesh grid
 ------------------------------------------------------------------------- */
int FixKineticsDiffusion::get_index(int grid) {
  int ix = floor(xGrid[grid][0] / stepx);
  int iy = floor(xGrid[grid][1] / stepy);
  int iz = floor(xGrid[grid][2] / stepz);

  int ind = iz * kinetics->nx * kinetics->ny + iy * kinetics->nx + ix;

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
    if (zbcflag == 0) {
      if (comm->nprocs < 2) { // periodic boundary conditions are already handled by cell communications for cases with more than 1 process 
	int zhiGrid = grid + nX * nY * nz;
	nuCell = nuPrev[zhiGrid];
      }
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
  else if (xGrid[grid][2] > zhi && !ghost[down]) {
    if (zbcflag == 0) {
      if (comm->nprocs < 2) {
	int zloGrid = grid - nX * nY * nz;
	nuCell = nuPrev[zloGrid];
      }
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
    if (ybcflag == 0) {
      if (comm->nprocs < 2) {
	int yhiGrid = grid + nX * ny;
	nuCell = nuPrev[yhiGrid];
      }
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
    if (ybcflag == 0) {
      if (comm->nprocs < 2) {
	int yloGrid = grid - nX * ny;
	nuCell = nuPrev[yloGrid];
      }
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
    if (xbcflag == 0) {
      if (comm->nprocs < 2) {
	int xhiGrid = grid + nx;
	nuCell = nuPrev[xhiGrid];
      }
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
    if (xbcflag == 0) {
      if (comm->nprocs < 2) {
	int xloGrid = grid - nx;
	nuCell = nuPrev[xloGrid];
      }
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

void FixKineticsDiffusion::compute_flux(double cellDNu, double &nuCell, double *nuPrev, double rateNu, int grid) {
  int lhs = grid - 1;   // x direction
  int rhs = grid + 1;  // x direction
  int bwd = grid - nX;  // y direction
  int fwd = grid + nX;  // y direction
  int down = grid - nX * nY; // z direction
  int up = grid + nX * nY;  // z dirction

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
    double uX = kinetics->fV[0][grid] * (nuPrev[rhs] - nuPrev[lhs]) / stepx;
    double uY = kinetics->fV[1][grid] * (nuPrev[fwd] - nuPrev[bwd]) / stepy;
    double uZ = kinetics->fV[2][grid] * (nuPrev[up] - nuPrev[down]) / stepz;

    res = (jX + jY + jZ + rateNu - uX - uY - uZ) * diffT;
  } else if (shearRate == 1) {
    int hgrid = grid;
    int deep = 0;

    while (!ghost[hgrid]) {
      hgrid = hgrid - nX * nY;
      deep++;
    }

    shear = shearRate * (deep * stepz - stepz / 2) * (nuPrev[rhs] - nuPrev[lhs]) / (2 * stepz);

    res = (jX + jY + jZ + rateNu - shear) * diffT;
  } else {
    res = (jX + jY + jZ + rateNu) * diffT;
  }

  //Updating the value: Ratesub*diffT + nuCell[cell](previous)
  nuCell = nuPrev[grid] + res;
}

/* ----------------------------------------------------------------------
 compare double values for equality
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

void FixKineticsDiffusion::update_nuS(){
  if ((xbcflag == 0 || xbcflag == 3) && (ybcflag == 0 || ybcflag == 3) && (zbcflag == 0 || zbcflag == 3)) {
    nuS = kinetics->nuS;
    nuR = kinetics->nuR;

    for (int nu = 1; nu < nnus+1; nu++) {
      for (int grid = 0; grid < nXYZ; grid++) {
        // transform nXYZ index to nuR index
        if (!ghost[grid]) {
          int ind = get_index(grid);
          double r = 0;

          if (unit == 0)
            r = nuR[nu][ind] * update->dt * kinetics->nevery * 1000;
          else
            r = nuR[nu][ind] * update->dt * kinetics->nevery;

          nuGrid[grid][nu] += r;

          if (nuGrid[grid][nu] < 0)
            nuGrid[grid][nu] = 0;

          if (unit == 0)
            nuS[nu][ind] = nuGrid[grid][nu] / 1000;
          else
            nuS[nu][ind] = nuGrid[grid][nu];
        }
      }
    }
  }
}
