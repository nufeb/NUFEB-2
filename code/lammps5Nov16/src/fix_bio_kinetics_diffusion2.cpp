/*
 * fix_kinetics/diffusion2.cpp
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
#include "fix_bio_kinetics_diffusion2.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "memory.h"
#include "update.h"
#include "variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion2::FixKineticsDiffusion2(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if (narg != 9) error->all(FLERR,"Not enough arguments in fix diffusion command");

  var = new char*[1];
  ivar = new int[1];

  if(strcmp(arg[3], "exp") == 0) sflag = 0;
  else if(strcmp(arg[3], "imp") == 0) sflag = 1;
  else error->all(FLERR,"Illegal PDE method command");

  for (int i = 0; i < 1; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  //set boundary condition flag:
  //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  if(strcmp(arg[5], "pp") == 0) xbcflag = 0;
  else if(strcmp(arg[5], "dd") == 0) xbcflag = 1;
  else if(strcmp(arg[5], "nd") == 0) xbcflag = 2;
  else if(strcmp(arg[5], "nn") == 0) xbcflag = 3;
  else if(strcmp(arg[5], "dn") == 0) xbcflag = 4;
  else error->all(FLERR,"Illegal x-axis boundary condition command");

  if(strcmp(arg[6], "pp") == 0) ybcflag = 0;
  else if(strcmp(arg[6], "dd") == 0) ybcflag = 1;
  else if(strcmp(arg[6], "nd") == 0) ybcflag = 2;
  else if(strcmp(arg[6], "nn") == 0) ybcflag = 3;
  else if(strcmp(arg[6], "dn") == 0) ybcflag = 4;
  else error->all(FLERR,"Illegal y-axis boundary condition command");

  if(strcmp(arg[7], "pp") == 0) zbcflag = 0;
  else if(strcmp(arg[7], "dd") == 0) zbcflag = 1;
  else if(strcmp(arg[7], "nd") == 0) zbcflag = 2;
  else if(strcmp(arg[7], "nn") == 0) zbcflag = 3;
  else if(strcmp(arg[7], "dn") == 0) zbcflag = 4;
  else error->all(FLERR,"Illegal z-axis boundary condition command");
  printf("%i, %i, %i \n", xbcflag, ybcflag, zbcflag);
  rflag = 0;
  rstep = atoi(arg[8]);
  if (rstep > 1) rflag = 1;

}

/* ---------------------------------------------------------------------- */

FixKineticsDiffusion2::~FixKineticsDiffusion2()
{
  int i;
  for (i = 0; i < 1; i++) {
    delete [] var[i];
  }

  delete [] var;
  delete [] ivar;

  memory->destroy(diffD);
  memory->destroy(xGrid);
  memory->destroy(nuGrid);
  memory->destroy(ghost);
}

/* ---------------------------------------------------------------------- */

int FixKineticsDiffusion2::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsDiffusion2::init()
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (!atom->radius_flag)
    error->all(FLERR,"Fix requires atom attribute diameter");

  for (int n = 0; n < 1; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix diffusion does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix diffusion is invalid style");
  }


  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required for running iBM simulation");

  bio = kinetics->bio;
  tol = input->variable->compute_equal(ivar[0]);

  //set diffusion grid size
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;
  nX = nx + 2;
  nY = ny + 2;
  nZ = nz + 2;

  nXYZ=(nX)*(nY)*(nZ);

  nnus = bio->nnus;
  iniS = bio->iniS;
  diffCoeff = bio->diffCoeff;

  diffD = memory->create(diffD,nnus+1,"diffusion:diffD");

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi-xlo)/nx;
  stepy = (yhi-ylo)/ny;
  stepz = (zhi-zlo)/nz;

  //if (!isEuqal(stepx, stepy, stepz)) error->all(FLERR,"Grid is not cubic");

  //inlet concentration, diffusion constant
  //and maximum boundary condition conc value
  for (int i = 1; i <= nnus; i++) {
    diffD[i] = diffCoeff[i];
  }

  xGrid =  memory->create(xGrid,nXYZ,3,"diffusion:xGrid");
  nuGrid = memory->create(nuGrid,nnus+1,nXYZ,"diffusion:nuGrid");
  ghost = memory->create(ghost,nXYZ,"diffusion:ghost");

  //initialise grids
  double i, j, k;
  int grid = 0;
  for (i = xlo - (stepx/2); i < xhi + stepx; i += stepx) {
    for (j = ylo - (stepy/2); j < yhi + stepy; j += stepy) {
      for (k = zlo - (stepz/2); k < zhi + stepz; k += stepz) {
        xGrid[grid][0] = i;
        xGrid[grid][1] = j;
        xGrid[grid][2] = k;
        //Initialise concentration values for ghost and std grids
        for (int nu = 1; nu <= nnus; nu++) {
          if (i < xlo) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][1] * 1000;
          } else if (i > xhi) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][2] * 1000;
          } else if (j < ylo) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][3] * 1000;
          } else if (j > yhi) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][4] * 1000;
          } else if (k < zlo) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][5] * 1000;
          } else if (k > zhi) {
            ghost[grid] = true;
            nuGrid[nu][grid] = iniS[nu][6] * 1000;
          } else {
            ghost[grid] = false;
            nuGrid[nu][grid] = iniS[nu][0] * 1000;
          }
        }
        grid++;
      }
    }
  }
}


/* ----------------------------------------------------------------------
  solve diffusion and reaction
------------------------------------------------------------------------- */

bool* FixKineticsDiffusion2::diffusion(bool *nuConv, int iter, double diffT, double **nuPrev)
{
  this->diffT = diffT;
  nuS = kinetics->nuS;
  nuR = kinetics->nuR;

  for (int i = 1; i <= nnus; i++) {
    if (bio->nuType[i] == 0 && diffCoeff[i] != 0 && !nuConv[i]) {
      double maxS = 0;

      // copy current concentrations
      for (int grid = 0; grid < nXYZ; grid++) {
        nuPrev[i][grid] = nuGrid[i][grid];
      }

      xbcm = iniS[i][1] * 1000;
      xbcp = iniS[i][2] * 1000;
      ybcm = iniS[i][3] * 1000;
      ybcp = iniS[i][4] * 1000;
      zbcm = iniS[i][5] * 1000;
      zbcp = iniS[i][6] * 1000;

      // solve diffusion and reaction
      for (int grid = 0; grid < nXYZ; grid++) {
        // transform nXYZ index to nuR index
        if (!ghost[grid]) {
          int ix = floor(xGrid[grid][0]/stepx);
          int iy = floor(xGrid[grid][1]/stepy);
          int iz = floor(xGrid[grid][2]/stepz);

          int ind = ix * ny * nz + iy * nz + iz;
         // printf("\n nuS1 = %e \n", nuS[i][ind]);
          compute_flux(diffD[i], nuGrid[i], nuPrev[i], nuR[i][ind], grid);
          printf("nuR[i][ind] = %e  \n", nuR[i][ind]);
          if (maxS < nuGrid[i][grid]) maxS = nuGrid[i][grid];

          if (nuGrid[i][grid] > 0) nuS[i][ind] = nuGrid[i][grid] / 1000;
          else {
            nuGrid[i][grid] = 0;
            nuS[i][ind] = 0;
          }
        } else compute_bc(nuGrid[i], nuPrev[i], grid);
      }

      bool conv = true;
      double r = -100;
      int gridx = 0;
      // check convergence criteria
      for (int grid = 0; grid < nXYZ; grid++) {
        if(!ghost[grid]){
          if (maxS == 0) maxS = 1;

          double rate = nuGrid[i][grid]/maxS;
          double prevRate = nuPrev[i][grid]/maxS;
          if (fabs(rate - prevRate) > r) r =fabs(rate - prevRate);
          if(fabs(rate - prevRate) >= tol)  {
            conv = false;
            gridx++;
            printf("nus = %e, prv = %e  \n", nuGrid[i][grid], nuPrev[i][grid]);
            //break;
          }
        }
      }
      nuConv[i] = conv;
      printf("grid = %i, %e  \n", gridx, r);
    }
  }
//
//  printf(" \n ");
//
//  for (int i = nZ-1; i>=0; i--){
//    printf(" \n ");
//    for (int j = 0; j<nY; j++){
//      printf("  ");
//      for (int k = 0; k<nX; k++){
//        int ind = k * nY * nZ + j * nZ + i;
//        //printf("%e ", nuGrid[1][ind]);
//        //printf("%d ", ghost[ind]);
//        printf("%d ", ind);
//      }
//    }
//  }
//  printf(" \n ");

  return nuConv;
}

/* ----------------------------------------------------------------------
  update boundary condition
------------------------------------------------------------------------- */

void FixKineticsDiffusion2::compute_bc(double *nuCell, double *nuPrev, int grid) {
  //for nx = ny = nz = 1 grids
  //2  11  20       5  14  23       8  17  26
  //1  10  19       4  13  22       7  16  25
  //0  9   18       3  12  21       6  15  24
  int lhs = grid - nZ*nY; // x direction
  int rhs = grid + nZ*nY; // x direction
  int bwd = grid - nZ; // y direction
  int fwd = grid + nZ; // y direction
  int down = grid - 1; // z direction
  int up = grid + 1; // z direction

  // assign values to the ghost-grids according to the boundary conditions.
  // If ghostcells are Neu then take the values equal from the adjacent cells.
  // if ghostcells are dirich then take the values equal to negative of the adjacent cells.
  // if ghostcells are mixed then zlo ghost cells are nuemann, zhi ghost cells are dirichlet, other four surfaces are periodic BC.
    // Low-z surface
  if (xGrid[grid][2] < zlo && !ghost[up]) {
    //0=PERIODIC-PERIODIC,  1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
    if (zbcflag == 0) {
      int zhiGrid = grid + nz;
      nuCell[grid] = nuPrev[zhiGrid];
    } else if (zbcflag == 1) {
      nuCell[grid] = 2*zbcm - nuPrev[up];
    } else if (zbcflag == 2) {
      nuCell[grid] = nuPrev[up];
    } else if (zbcflag == 3) {
      nuCell[grid] = nuPrev[up];
    } else if (zbcflag == 4) {
      nuCell[grid] = 2*zbcm - nuPrev[up];
    }
  }
  // high-z surface
  else if (xGrid[grid][2] > zhi && !ghost[down]) {
    if (zbcflag == 0) {
      int zloGrid = grid - nz;
      nuCell[grid] = nuPrev[zloGrid];
    } else if (zbcflag == 1) {
      nuCell[grid] = 2*zbcp - nuPrev[down];
    } else if (zbcflag == 2) {
      nuCell[grid] = 2*zbcp - nuPrev[down];
    } else if (zbcflag == 3) {
      nuCell[grid] = nuPrev[down];
    } else if (zbcflag == 4) {
      nuCell[grid] = nuPrev[down];
    }
  }
  // low-y surface
  else if (xGrid[grid][1] < ylo && !ghost[fwd]) {
    if (ybcflag == 0) {
      int yhiGrid = grid + nZ*ny;
      nuCell[grid] = nuPrev[yhiGrid];
    } else if (ybcflag == 1) {
      nuCell[grid] = 2*ybcm - nuPrev[fwd];
    } else if (ybcflag == 2) {
      nuCell[grid] = nuPrev[fwd];
    } else if (ybcflag == 3) {
      nuCell[grid] = nuPrev[fwd];
    } else if (ybcflag == 4) {
      nuCell[grid] = 2*ybcm - nuPrev[fwd];
    }
  }
  // high-y surface
  else if (xGrid[grid][1] > yhi && !ghost[bwd]) {
    if (ybcflag == 0) {
      int yloGrid = grid - nZ*ny;
      nuCell[grid] = nuPrev[yloGrid];
    } else if (ybcflag == 1) {
      nuCell[grid] = 2*ybcp - nuPrev[bwd];
    } else if (ybcflag == 2) {
      nuCell[grid] = 2*ybcp - nuPrev[bwd];
    } else if (ybcflag == 3) {
      nuCell[grid] = nuPrev[bwd];
    } else if (ybcflag == 4) {
      nuCell[grid] = nuPrev[bwd];
    }
  }
  // low-x surface
  else if (xGrid[grid][0] < xlo && !ghost[rhs]) {
    if (xbcflag == 0) {
      int xhiGrid = grid + nY*nZ*nx;
      nuCell[grid] = nuPrev[xhiGrid];
    } else if (xbcflag == 1) {
      nuCell[grid] = 2*xbcm - nuPrev[rhs];
    } else if (xbcflag == 2) {
      nuCell[grid] = nuPrev[rhs];
    } else if (xbcflag == 3) {
      nuCell[grid] = nuPrev[rhs];
    } else if (xbcflag == 4) {
      nuCell[grid] = 2*xbcm - nuPrev[rhs];
    }
  }
  // high-x surface
  else if (xGrid[grid][0] > xhi && !ghost[lhs]) {
    if (xbcflag == 0) {
      int xloGrid = grid - nY*nZ*nx;
      nuCell[grid] = nuPrev[xloGrid];
    } else if (xbcflag == 1) {
      nuCell[grid] = 2*xbcp - nuPrev[lhs];
    } else if (xbcflag == 2) {
      nuCell[grid] = 2*xbcp - nuPrev[lhs];
    } else if (xbcflag == 3) {
      nuCell[grid] = nuPrev[lhs];
    } else if (xbcflag == 4) {
      nuCell[grid] = nuPrev[lhs];
    }
  }
}

/* ----------------------------------------------------------------------
  update non-ghost grids
------------------------------------------------------------------------- */

void FixKineticsDiffusion2::compute_flux(double cellDNu, double *nuCell, double *nuPrev, double rateNu, int grid) {
  //for nx = ny = nz = 1 grids
  //2  11  20       5  14  23       8  17  26
  //1  10  19       4  13  22       7  16  25
  //0  9   18       3  12  21       6  15  24
  int lhs = grid - nZ*nY; // x direction
  int rhs = grid + nZ*nY; // x direction
  int bwd = grid - nZ; // y direction
  int fwd = grid + nZ; // y direction
  int down = grid - 1; // z direction
  int up = grid + 1; // z direction

  double jRight = cellDNu*(nuPrev[rhs] - nuPrev[grid])/stepx;
  double jLeft = cellDNu*(nuPrev[grid] - nuPrev[lhs])/stepx;
  double jX = (jRight - jLeft)/stepx;

  double jForward = cellDNu*(nuPrev[fwd] - nuPrev[grid])/stepy;
  double jBackward = cellDNu*(nuPrev[grid] - nuPrev[bwd])/stepy;
  double jY = (jForward - jBackward)/stepy;

  double jUp = cellDNu*(nuPrev[up] - nuPrev[grid])/stepz;
  double jDown = cellDNu*(nuPrev[grid] - nuPrev[down])/stepz;
  double jZ = (jUp - jDown)/stepz;

  // Adding fluxes in all the directions and the uptake rate (RHS side of the equation)
  double Ratesub = jX + jY + jZ + rateNu;
  //Updating the value: Ratesub*diffT + nuCell[cell](previous)
  nuCell[grid] = nuPrev[grid] + Ratesub*diffT;
//  printf("before = %e \n", nuPrev[grid]);
//  printf("after = %e \n", nuCell[grid]);
}


/* ----------------------------------------------------------------------
  compare double values for equality
------------------------------------------------------------------------- */

bool FixKineticsDiffusion2::isEuqal(double a, double b, double c)
{
  double epsilon = 1e-10;
  if ((fabs(a - b) > epsilon)|| (fabs(a - b) > epsilon) || (fabs(a - b) > epsilon))
    return false;

  return true;
}
