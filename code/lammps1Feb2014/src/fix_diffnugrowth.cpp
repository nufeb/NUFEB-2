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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_diffnugrowth.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDiffNuGrowth::FixDiffNuGrowth(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 48) error->all(FLERR,"Not enough arguments in fix diff growth command");

  nevery = force->inumeric(FLERR,arg[3]);
  diffevery = force->inumeric(FLERR,arg[4]);
  outputevery = force->inumeric(FLERR,arg[47]);
  if (nevery < 0 || diffevery < 0 || outputevery < 0) error->all(FLERR,"Illegal fix growth command");

  var = new char*[38];
  ivar = new int[38];

  for (int i = 0; i < 33; i++) {
    int n = strlen(&arg[5+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[5+i][2]);
  }

  //BC concentration
	for(int i = 33; i < 38; i++){
		 int n = strlen(&arg[9+i][2]) + 1;
		 var[i] = new char[n];
		 strcpy(var[i],&arg[9+i][2]);
	}

  if(strcmp(arg[41], "dirich") == 0) bflag = 1;
  else if(strcmp(arg[41], "neu") == 0) bflag = 2;
  else if(strcmp(arg[41], "mixed") == 0) bflag = 3;
  else error->all(FLERR,"Illegal boundary condition command");

  nx = atoi(arg[38]);
  ny = atoi(arg[39]);
  nz = atoi(arg[40]);

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
}

/* ---------------------------------------------------------------------- */

FixDiffNuGrowth::~FixDiffNuGrowth()
{
  int i;
  for (i = 0; i < 38; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
  delete [] xCell;
  delete [] yCell;
  delete [] zCell;
  delete [] cellVol;
  delete [] ghost;
  delete [] subCell;
  delete [] o2Cell;
  delete [] nh4Cell;
  delete [] no2Cell;
  delete [] no3Cell;
}

/* ---------------------------------------------------------------------- */

int FixDiffNuGrowth::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDiffNuGrowth::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix growth requires atom attribute diameter");

  for (int n = 0; n < 38; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix nugrowth does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix nugrowth is invalid style");
  }

  //initial concentrations of non-ghost cells
  initsub = input->variable->compute_equal(ivar[28]);
  inito2 = input->variable->compute_equal(ivar[29]);
  initno2 = input->variable->compute_equal(ivar[30]);
  initno3 = input->variable->compute_equal(ivar[31]);
  initnh4 = input->variable->compute_equal(ivar[32]);

  //initial concentrations of boundary
	subBC = input->variable->compute_equal(ivar[33]);
	o2BC = input->variable->compute_equal(ivar[34]);
	no2BC = input->variable->compute_equal(ivar[35]);
	no3BC = input->variable->compute_equal(ivar[36]);
	nh4BC = input->variable->compute_equal(ivar[37]);

  //total numbers of cells (ghost + non-ghost)
  numCells = (nx+2)*(ny+2)*(nz+2);

  xCell = new double[numCells];
  yCell = new double[numCells];
  zCell = new double[numCells];
  cellVol = new double[numCells];
  ghost = new bool[numCells];
  subCell = new double[numCells];
  o2Cell = new double[numCells];
  nh4Cell = new double[numCells];
  no2Cell = new double[numCells];
  no3Cell = new double[numCells];

  xstep = (xhi - xlo) / nx;
  ystep = (yhi - ylo) / ny;
  zstep = (zhi - zlo) / nz;

  //initialise cells
  double i, j, k;
  int cell = 0;
  for (i = xlo - (xstep/2); i < xhi + xstep; i += xstep) {
    for (j = ylo - (ystep/2); j < yhi + ystep; j += ystep) {
      for (k = zlo - (zstep/2); k < zhi + zstep; k += zstep) {
        xCell[cell] = i;
        yCell[cell] = j;
        zCell[cell] = k;
        cellVol[cell] = xstep * ystep * zstep;
        ghost[cell] = false;
        //Initialise concentration values for ghost and std cells
        if (i < xlo || i > xhi || j < ylo ||
        	j > yhi || k < zlo || k > zhi) {
        		ghost[cell] = true;

						subCell[cell] = subBC;
						o2Cell[cell] = o2BC;
						no2Cell[cell] = no2BC;
						no3Cell[cell] = no3BC;
						nh4Cell[cell] = nh4BC;
        }else{
            subCell[cell] = initsub;
            o2Cell[cell] = inito2;
            no2Cell[cell] = initno2;
            no3Cell[cell] = initno3;
            nh4Cell[cell] = initnh4;
        }
        cell++;
      }
    }
  }

  //create folders for concentration data
  int status1, status2;

  status1 = system("rm -rf sub o2 no2 no3 nh4");
	status2 = system("mkdir sub o2 no2 no3 nh4");

	if(status1 < 0 || status2 < 0){
		 error->all(FLERR,"Fail to create concentration output dir");
	}
}

/* ---------------------------------------------------------------------- */

void FixDiffNuGrowth::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  change_dia();
}

/* ---------------------------------------------------------------------- */

void FixDiffNuGrowth::change_dia()
{
  modify->clearstep_compute();

  double KsHET = input->variable->compute_equal(ivar[0]);
  double Ko2HET = input->variable->compute_equal(ivar[1]);
  double Kno2HET = input->variable->compute_equal(ivar[2]);
  double Kno3HET = input->variable->compute_equal(ivar[3]);
  double Knh4AOB = input->variable->compute_equal(ivar[4]);
  double Ko2AOB = input->variable->compute_equal(ivar[5]);
  double Kno2NOB = input->variable->compute_equal(ivar[6]);
  double Ko2NOB = input->variable->compute_equal(ivar[7]);
  double MumHET = input->variable->compute_equal(ivar[8]);
  double MumAOB = input->variable->compute_equal(ivar[9]);
  double MumNOB = input->variable->compute_equal(ivar[10]);
  double etaHET = input->variable->compute_equal(ivar[11]);
  double bHET = input->variable->compute_equal(ivar[12]); // R6
  double bAOB = input->variable->compute_equal(ivar[13]); // R7
  double bNOB = input->variable->compute_equal(ivar[14]); // R8
  double bEPS = input->variable->compute_equal(ivar[15]); // R9
  double YHET = input->variable->compute_equal(ivar[16]);
  double YAOB = input->variable->compute_equal(ivar[17]);
  double YNOB = input->variable->compute_equal(ivar[18]);
  double YEPS = input->variable->compute_equal(ivar[19]);
  double Y1 = input->variable->compute_equal(ivar[20]);
  double EPSdens = input->variable->compute_equal(ivar[21]);
  double Do2 = input->variable->compute_equal(ivar[22]);
  double Dnh4 = input->variable->compute_equal(ivar[23]);
  double Dno2 = input->variable->compute_equal(ivar[24]);
  double Dno3 = input->variable->compute_equal(ivar[25]);
  double Ds = input->variable->compute_equal(ivar[26]);
  double diffT = input->variable->compute_equal(ivar[27]);

  double density;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outerMass = atom->outerMass;
  double *outerRadius = atom->outerRadius;
  double *sub = atom->sub;
  double *o2 = atom->o2;
  double *nh4 = atom->nh4;
  double *no2 = atom->no2;
  double *no3 = atom->no3;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

  int* cellIn = new int[nall];
  double* xHET = new double[numCells];
  double* xAOB = new double[numCells];
  double* xNOB = new double[numCells];
  double* xEPS = new double[numCells];
  double* xTot = new double[numCells];

  for (int cell = 0; cell < numCells; cell++) {
  	xHET[cell] = 0.0;
  	xAOB[cell] = 0.0;
  	xNOB[cell] = 0.0;
  	xEPS[cell] = 0.0;
  	xTot[cell] = 0.0;
  }

  double* R1 = new double[numCells];
  double* R2 = new double[numCells];
  double* R3 = new double[numCells];
  double* R4 = new double[numCells];
  double* R5 = new double[numCells];
  // double R6[numCells] = bHET;
  // double R7[numCells] = bAOB;
  // double R8[numCells] = bNOB;
  // double R9[numCells] = bEPS;
  double* Rs = new double[numCells];
  double* Ro2 = new double[numCells];
  double* Rnh4 = new double[numCells];
  double* Rno2 = new double[numCells];
  double* Rno3 = new double[numCells];

  double* cellDo2 = new double[numCells];
  double* cellDnh4 = new double[numCells];
  double* cellDno2 = new double[numCells];
  double* cellDno3  = new double[numCells];
  double* cellDs  = new double[numCells];

  int grid = 0;
  // Figure out which cell each particle is in
  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      double gHET = 0;
      double gAOB = 0;
      double gNOB = 0;
      double gEPS = 0;
      if (type[i] == 1) {
        gHET = 1;
       // fprintf(stdout, "mass = %e\n", atom->rmass[i]);
      }
      if (type[i] == 2) {
        gAOB = 1;
      }
      if (type[i] == 3) {
        gNOB = 1;
      }
      if (type[i] == 4) {
        gEPS = 1;
      }

      bool allocate = false;
      for (int j = 0; j < numCells; j ++) {
        if ((xCell[j] - xstep/2) <= atom->x[i][0] &&
            (xCell[j] + xstep/2) >= atom->x[i][0] &&
            (yCell[j] - ystep/2) <= atom->x[i][1] &&
            (yCell[j] + ystep/2) >= atom->x[i][1] &&
            (zCell[j] - zstep/2) <= atom->x[i][2] &&
            (zCell[j] + zstep/2) >= atom->x[i][2]) {
          cellIn[i] = j;
          xHET[j] += (gHET * rmass[i])/cellVol[j];
          xAOB[j] += (gAOB * rmass[i])/cellVol[j];
        //  xNOB[j] += rmass[i]/cellVol[j];
          xNOB[j] += (gNOB *rmass[i])/cellVol[j];
          xEPS[j] += (gEPS * rmass[i])/cellVol[j];
          xTot[j] += rmass[i]/cellVol[j];
          allocate = true;
          break;
        }
      }
      if(!allocate)
      	error->all(FLERR,"Fail to allocate grid.");
    }
  }

  //initialize values
  for (int cell = 0; cell < numCells; cell++) {
    R1[cell] = MumHET*(subCell[cell]/(KsHET+subCell[cell]))*(o2Cell[cell]/(Ko2HET+o2Cell[cell]));
    R2[cell] = MumAOB*(nh4Cell[cell]/(Knh4AOB+nh4Cell[cell]))*(o2Cell[cell]/(Ko2AOB+o2Cell[cell]));
    R3[cell] = MumNOB*(no2Cell[cell]/(Kno2NOB+no2Cell[cell]))*(o2Cell[cell]/(Ko2NOB+o2Cell[cell]));
    R4[cell] = etaHET*MumHET*(subCell[cell]/(KsHET+subCell[cell]))*(no3Cell[cell]/(Kno3HET+no3Cell[cell]))*(Ko2HET/(Ko2HET+o2Cell[cell]));
    R5[cell] = etaHET*MumHET*(subCell[cell]/(KsHET+subCell[cell]))*(no2Cell[cell]/(Kno2HET+no2Cell[cell]))*(Ko2HET/(Ko2HET+o2Cell[cell]));

    if(!(update->ntimestep % diffevery)){
			Rs[cell] = ( (-1/YHET) * ( (R1[cell]+R4[cell]+R5[cell]) * xHET[cell] ) ) + ( (1-Y1) * ( bHET*xHET[cell]+bAOB*xAOB[cell]+bNOB*xNOB[cell] ) ) + ( bEPS*xEPS[cell] );
			Ro2[cell] = (-((1-YHET-YEPS)/YHET)*R1[cell]*xHET[cell])-(((3.42-YAOB)/YAOB)*R2[cell]*xAOB[cell])-(((1.15-YNOB)/YNOB)*R3[cell]*xNOB[cell]);
			Rnh4[cell] = -(1/YAOB)*R2[cell]*xAOB[cell];
			Rno2[cell] = ((1/YAOB)*R2[cell]*xAOB[cell])-((1/YNOB)*R3[cell]*xNOB[cell])-(((1-YHET-YEPS)/(1.17*YHET))*R5[cell]*xHET[cell]);
			Rno3[cell] = ((1/YNOB)*R3[cell]*xNOB[cell])-(((1-YHET-YEPS)/(2.86*YHET))*R4[cell]*xHET[cell]);

			double diffusionFunction = 1 - ((0.43 * pow(xTot[cell], 0.92))/(11.19+0.27*pow(xTot[cell], 0.99)));

			cellDo2[cell] = diffusionFunction * Do2;
			cellDnh4[cell] = diffusionFunction * Dnh4;
			cellDno2[cell] = diffusionFunction * Dno2;
			cellDno3[cell] = diffusionFunction * Dno3;
			cellDs[cell] = diffusionFunction * Ds;
    }
  }
  if(!(update->ntimestep % diffevery)) {

  	double* subPrev  = new double[numCells];
  	double* o2Prev  = new double[numCells];
  	double* no2Prev  = new double[numCells];
  	double* no3Prev  = new double[numCells];
  	double* nh4Prev  = new double[numCells];

		bool subConvergence = false;
		bool o2Convergence = false;
		bool no2Convergence = false;
		bool no3Convergence = false;
		bool nh4Convergence = false;

		bool convergence = false;

		int iteration = 0;

		double tol = 1e-6; // Tolerance for convergence criteria for nutrient balance equation

		// Outermost while loop for the convergence criterion
		while (!convergence) {
			iteration ++;

			for (int cell = 0; cell < numCells; cell++) {
				subPrev[cell] = subCell[cell];
				o2Prev[cell] = o2Cell[cell];
				nh4Prev[cell] = nh4Cell[cell];
				no2Prev[cell] = no2Cell[cell];
				no3Prev[cell] = no3Cell[cell];
			}

			for (int cell = 0; cell < numCells; cell++) {
				R1[cell] = MumHET*(subPrev[cell]/(KsHET+subPrev[cell]))*(o2Prev[cell]/(Ko2HET+o2Prev[cell]));
				R2[cell] = MumAOB*(nh4Prev[cell]/(Knh4AOB+nh4Prev[cell]))*(o2Prev[cell]/(Ko2AOB+o2Prev[cell]));
				R3[cell] = MumNOB*(no2Prev[cell]/(Kno2NOB+no2Prev[cell]))*(o2Prev[cell]/(Ko2NOB+o2Prev[cell]));
				R4[cell] = etaHET*MumHET*(subPrev[cell]/(KsHET+subPrev[cell]))*(no3Prev[cell]/(Kno3HET+no3Prev[cell]))*(Ko2HET/(Ko2HET+o2Prev[cell]));
				R5[cell] = etaHET*MumHET*(subPrev[cell]/(KsHET+subPrev[cell]))*(no2Prev[cell]/(Kno2HET+no2Prev[cell]))*(Ko2HET/(Ko2HET+o2Prev[cell]));

				Rs[cell] = ((-1/YHET) * ( (R1[cell]+R4[cell]+R5[cell]) * xHET[cell] ) ) + ( (1-Y1) * ( bHET*xHET[cell]+bAOB*xAOB[cell]+bNOB*xNOB[cell] ) ) + ( bEPS*xEPS[cell] );
				Ro2[cell] = (-((1-YHET-YEPS)/YHET)*R1[cell]*xHET[cell])-(((3.42-YAOB)/YAOB)*R2[cell]*xAOB[cell])-(((1.15-YNOB)/YNOB)*R3[cell]*xNOB[cell]);
				Rnh4[cell] = -(1/YAOB)*R2[cell]*xAOB[cell];
				Rno2[cell] = ((1/YAOB)*R2[cell]*xAOB[cell])-((1/YNOB)*R3[cell]*xNOB[cell])-(((1-YHET-YEPS)/(1.17*YHET))*R5[cell]*xHET[cell]);
				Rno3[cell] = ((1/YNOB)*R3[cell]*xNOB[cell])-(((1-YHET-YEPS)/(2.86*YHET))*R4[cell]*xHET[cell]);

	    	if(!subConvergence) computeFlux(cellDs, subCell, subPrev, subBC, Rs[cell], diffT, cell);
				if(!o2Convergence) computeFlux(cellDo2, o2Cell, o2Prev, o2BC, Ro2[cell], diffT, cell);
				if(!nh4Convergence) computeFlux(cellDnh4, nh4Cell, nh4Prev, nh4BC, Rnh4[cell], diffT, cell);
				if(!no2Convergence) computeFlux(cellDno2, no2Cell, no2Prev, no2BC, Rno2[cell], diffT, cell);
				if(!no3Convergence) computeFlux(cellDno3, no3Cell, no3Prev, no3BC, Rno3[cell], diffT, cell);
				}


			if(isConvergence(subCell, subPrev, subBC, tol))	subConvergence = true;
			if(isConvergence(o2Cell, o2Prev, o2BC, tol))	o2Convergence = true;
			if(isConvergence(nh4Cell, nh4Prev, nh4BC, tol)) nh4Convergence = true;
			if(isConvergence(no2Cell, no2Prev, no2BC, tol)) no2Convergence = true;
			if(isConvergence(no3Cell, no3Prev, no3BC, tol)) no3Convergence = true;

			if((subConvergence && o2Convergence && nh4Convergence && no2Convergence && no3Convergence) || iteration == 10000) {
				convergence = true;
			}
		}

		fprintf(stdout, "Number of iterations:  %i\n", iteration);

	  delete [] subPrev;
	  delete [] o2Prev;
	  delete [] nh4Prev;
	  delete [] no2Prev;
	  delete [] no3Prev;
  }

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      double gHET = 0;
      double gAOB = 0;
      double gNOB = 0;
      double gEPS = 0;

      if (type[i] == 1) {
        gHET = 1;
      }
      if (type[i] == 2) {
        gAOB = 1;
      }
      if (type[i] == 3) {
        gNOB = 1;
      }
      if (type[i] == 4) {
        gEPS = 1;
      }

      sub[i] = subCell[cellIn[i]];
      o2[i] = o2Cell[cellIn[i]];
      nh4[i] = nh4Cell[cellIn[i]];
      no2[i] = no2Cell[cellIn[i]];
      no3[i] = no3Cell[cellIn[i]];

      double value = update->dt * (gHET*(R1[cellIn[i]]+R4[cellIn[i]]+R5[cellIn[i]]) + gAOB*R2[cellIn[i]] + gNOB*R3[cellIn[i]] - gEPS*bEPS);

      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      double oldMass = rmass[i];
      rmass[i] = rmass[i]*(1 + (value*nevery));
      //if set critical minimal value to 0, may get double out of range error.
      if (rmass[i] <= 1e-20) {
        rmass[i] = oldMass;
      }
      
      double value2 = update->dt * (YEPS/YHET)*(R1[cellIn[i]]+R4[cellIn[i]]+R5[cellIn[i]]);
      double oldRadius = radius[i];
      if (type[i] == 1) {
        outerMass[i] = (((4.0*MY_PI/3.0)*((outerRadius[i]*outerRadius[i]*outerRadius[i])-(radius[i]*radius[i]*radius[i])))*EPSdens)+(value2*nevery*rmass[i]);

        outerRadius[i] = pow((3.0/(4.0*MY_PI))*((rmass[i]/density)+(outerMass[i]/EPSdens)),(1.0/3.0));
        radius[i] = pow((3.0/(4.0*MY_PI))*(rmass[i]/density),(1.0/3.0));
      }
      else {
        radius[i] = pow((3.0/(4.0*MY_PI))*(rmass[i]/density),(1.0/3.0));
        outerMass[i] = 0.0;
        outerRadius[i] = radius[i];
      }
    }
  }

	//output concentration
	if(!(update->ntimestep % outputevery)){
	  outputData(outputevery,1);
	  outputData(outputevery,2);
	  outputData(outputevery,3);
	  outputData(outputevery,4);
	  outputData(outputevery,5);
	}

  modify->addstep_compute(update->ntimestep + nevery);

  delete [] cellIn;
  delete [] xHET;
  delete [] xAOB;
  delete [] xNOB;
  delete [] xEPS;
  delete [] xTot;

  delete [] R1;
  delete [] R2;
  delete [] R3;
  delete [] R4;
  delete [] R5;

  delete [] Rs;
  delete [] Ro2;
  delete [] Rnh4;
  delete [] Rno3;
  delete [] Rno2;

  delete [] cellDo2;
  delete [] cellDnh4;
  delete [] cellDno2;
  delete [] cellDno3;
  delete [] cellDs;
}

bool FixDiffNuGrowth::isConvergence(double *nuCell, double *prevNuCell, double nuBC, double tol) {
	for(int cell = 0; cell < numCells; cell++){
		if(!ghost[cell]){

			double rate = nuCell[cell]/nuBC;
			double prevRate = prevNuCell[cell]/nuBC;

			if(fabs(rate - prevRate) >= tol) return false;
		}
	}
	return true;
}

void FixDiffNuGrowth::computeFlux(double *cellDNu, double *nuCell, double *nuPrev, double nuBC, double rateNu, double diffT, int cell) {
 // fprintf(stdout, "nuCell= %f, nuPrev = %f, rateNu = %f\n", nuCell[cell], nuPrev[cell], rateNu);
	//for nx = ny = nz = 1 grids
	//2  11  20 			5  14  23				8  17  26
	//1  10  19       4  13  22       7  16  25
	//0  9   18       3  12  21       6  15  24
	int leftCell = cell - (nz+2)*(ny+2); // x direction
	int rightCell = cell + (nz+2)*(ny+2); // x direction
	int downCell = cell - (nz+2); // y direction
	int upCell = cell + (nz+2); // y direction
	int backwardCell = cell - 1; // z direction
	int forwardCell = cell + 1; // z direction

	// assign values to the ghost-cells according to the boundary conditions.
	// If ghostcells are Neu then take the values equal from the adjacent cells.
	// if ghostcells are dirich then take the values equal to negative of the adjacent cells.
	// if ghostcells are mixed then zlo ghost cells are nuemann, zhi ghost cells are dirichlet, other four surfaces are periodic BC.
	if (ghost[cell]) {
		// fprintf(stdout, "Ghost Cell: %i\n", cell);
		if (zCell[cell] < zlo && !ghost[forwardCell]) {
			if (bflag == 1) {
				nuCell[cell] = 2*nuBC - nuPrev[forwardCell];
			} else if(bflag == 2 || bflag == 3) {
				nuCell[cell] = nuPrev[forwardCell];
			}
		}
		else if (zCell[cell] > zhi && !ghost[backwardCell]) {
			if (bflag == 1 || bflag == 3) {
				nuCell[cell] = 2*nuBC - nuPrev[backwardCell];
			} else if (bflag == 2) {
				nuCell[cell] = nuPrev[backwardCell];
			}
		}
		else if (yCell[cell] < ylo && !ghost[upCell]) {
			if (bflag == 1) {
				nuCell[cell] = 2*nuBC - nuPrev[upCell];
			} else if (bflag == 2) {
				nuCell[cell] = nuPrev[upCell];
			} else if (bflag == 3) {
				int yhiCell = cell + (nz+2)*ny;
				nuCell[cell] = nuPrev[yhiCell];
			}
		}
		else if (yCell[cell] > yhi && !ghost[downCell]) {
			if (bflag == 1) {
				nuCell[cell] = 2*nuBC - nuPrev[downCell];
			} else if (bflag == 2) {
				nuCell[cell] = nuPrev[downCell];
			} else if (bflag == 3) {
				int yloCell = cell - (nz+2)*ny;
				nuCell[cell] = nuPrev[yloCell];
			}
		}
		else if (xCell[cell] < xlo && !ghost[rightCell]) {
			if (bflag == 1) {
				nuCell[cell] = 2*nuBC - nuPrev[rightCell];
			} else if (bflag == 2) {
				nuCell[cell] = nuPrev[rightCell];
			} else if (bflag == 3) {
				int xhiCell = cell + (ny+2)*(nz+2)*nx;
				nuCell[cell] = nuPrev[xhiCell];
			}
		}
		else if (xCell[cell] > xhi && !ghost[leftCell]) {
			if (bflag == 1) {
				nuCell[cell] = 2*nuBC - nuPrev[leftCell];
			} else if (bflag == 2) {
				nuCell[cell] = nuPrev[leftCell];
			} else if (bflag == 3) {
				int xloCell = cell - (ny+2)*(nz+2)*nx;
				nuCell[cell] = nuPrev[xloCell];
			}
		}
	}
	else {
		double dRight = (cellDNu[cell] + cellDNu[rightCell]) / 2;
		double jRight = dRight*(nuPrev[rightCell] - nuPrev[cell])/xstep;
		double dLeft = (cellDNu[cell] + cellDNu[leftCell]) / 2;
		double jLeft = dLeft*(nuPrev[cell] - nuPrev[leftCell])/xstep;
		double jX = (jRight - jLeft)/xstep;

		double dUp = (cellDNu[cell] + cellDNu[upCell]) / 2;
		double jUp = dUp*(nuPrev[upCell] - nuPrev[cell])/ystep;
		double dDown = (cellDNu[cell] + cellDNu[downCell]) / 2;
		double jDown = dDown*(nuPrev[cell] - nuPrev[downCell])/ystep;
		double jY = (jUp - jDown)/ystep;

		double dForward = (cellDNu[cell] + cellDNu[forwardCell]) / 2;
		double jForward = dForward*(nuPrev[forwardCell] - nuPrev[cell])/zstep;
		double dBackward = (cellDNu[cell] + cellDNu[backwardCell]) / 2;
		double jBackward = dBackward*(nuPrev[cell] - nuPrev[backwardCell])/zstep;
		double jZ = (jForward - jBackward)/zstep;

		// Adding fluxes in all the directions and the uptake rate (RHS side of the equation)
		double Ratesub = jX + jY + jZ + rateNu;
		//Updating the value: Ratesub*diffT + nuCell[cell](previous)
		nuCell[cell] += Ratesub*diffT;

		if(nuCell[cell] <= 1e-20){
			nuCell[cell] = 1e-20;
		}
	}
}


void FixDiffNuGrowth::outputData(int every, int n){
  if(!(update->ntimestep%every)){
	  FILE* pFile;
	  std::string str;
	  std::ostringstream stm;
	  stm << update->ntimestep;

	  switch(n) {
	  case 1 :
	  	str = "./sub/sub.csv."+ stm.str();
	  	break;
	  case 2 :
	  	str = "./o2/o2.csv."+ stm.str();
	  	break;
	  case 3 :
	  	str = "./no2/no2.csv."+ stm.str();
	  	break;
	  case 4 :
	  	str = "./no3/no3.csv."+ stm.str();
	  	break;
	  case 5 :
	  	str = "./nh4/nh4.csv."+ stm.str();
	  	break;
	  }

	  pFile = fopen (str.c_str(), "w");

	  fprintf(pFile, ",x,y,z,scalar,1,1,1,0.5\n");
	  for(int i = 0; i < numCells; i++){
		  if(!ghost[i]){

			  switch(n) {
			  case 1 :
			  	fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, xCell[i], yCell[i], zCell[i], subCell[i]);
			  	break;
			  case 2 :
			  	fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, xCell[i], yCell[i], zCell[i], o2Cell[i]);
			  	break;
			  case 3 :
			  	fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, xCell[i], yCell[i], zCell[i], no2Cell[i]);
			  	break;
			  case 4 :
			  	fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, xCell[i], yCell[i], zCell[i], no3Cell[i]);
			  	break;
			  case 5 :
			  	fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, xCell[i], yCell[i], zCell[i], nh4Cell[i]);
			  	break;
			  }
		  }
	  }

	  fclose(pFile);
  }
}
