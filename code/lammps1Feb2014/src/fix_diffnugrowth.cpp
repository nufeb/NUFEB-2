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

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDiffNuGrowth::FixDiffNuGrowth(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 48) error->all(FLERR,"Not enough arguments in fix diff growth command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix growth command");

  var = new char*[33];
  ivar = new int[33];

  int i;
  for (i = 0; i < 33; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  nx = atoi(arg[37]);
  ny = atoi(arg[38]);
  nz = atoi(arg[39]);

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

  xloDirch = false;
  xhiDirch = false;
  yloDirch = false;
  yhiDirch = false;
  zloDirch = false;
  zhiDirch = false;

  if (strcmp(arg[40],"dirch") == 0) {
  	i = 41;
  	while (strcmp(arg[i],"neu") != 0) {
  		if (strcmp(arg[i],"xlo") == 0) {
  			xloDirch = true;
  		}
  		else if (strcmp(arg[i],"xhi") == 0) {
  			xhiDirch = true;
  		}
  		else if (strcmp(arg[i],"ylo") == 0) {
  			yloDirch = true;
  		}
  		else if (strcmp(arg[i],"yhi") == 0) {
  			yhiDirch = true;
  		}
  		else if (strcmp(arg[i],"zlo") == 0) {
  			zloDirch = true;
  		}
  		else if (strcmp(arg[i],"zhi") == 0) {
  			zhiDirch = true;
  		}
  		i++;
  	}
  }
  else if (strcmp(arg[40],"neu") == 0) {
  	i = 41;
  	while (strcmp(arg[i],"dirch") != 0) {
  		i++;
  	}
  	for (;i < narg; i++) {
  		if (strcmp(arg[i],"xlo") == 0) {
  			xloDirch = true;
  		}
  		else if (strcmp(arg[i],"xhi") == 0) {
  			xhiDirch = true;
  		}
  		else if (strcmp(arg[i],"ylo") == 0) {
  			yloDirch = true;
  		}
  		else if (strcmp(arg[i],"yhi") == 0) {
  			yhiDirch = true;
  		}
  		else if (strcmp(arg[i],"zlo") == 0) {
  			zloDirch = true;
  		}
  		else if (strcmp(arg[i],"zhi") == 0) {
  			zhiDirch = true;
  		}
  	}
  }
}

/* ---------------------------------------------------------------------- */

FixDiffNuGrowth::~FixDiffNuGrowth()
{
  int i;
  for (i = 0; i < 33; i++) {
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

  int n;
  for (n = 0; n < 33; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix nugrowth does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix nugrowth is invalid style");
  }

  double sub = input->variable->compute_equal(ivar[28]);
  double o2 = input->variable->compute_equal(ivar[29]);
  double no2 = input->variable->compute_equal(ivar[30]);
  double no3 = input->variable->compute_equal(ivar[31]);
  double nh4 = input->variable->compute_equal(ivar[32]);

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

  double i, j, k;
  int cell = 0;
  for (i = xlo - (xstep/2); i < xhi + xstep; i += xstep) {
    for (j = ylo - (ystep/2); j < yhi + ystep; j += ystep) {
      for (k = zlo - (zstep/2); k < zhi + zstep; k += zstep) {
        xCell[cell] = i;
        yCell[cell] = j;
        zCell[cell] = k;
        subCell[cell] = sub;
        o2Cell[cell] = o2;
        no2Cell[cell] = no2;
        no3Cell[cell] = no3;
        nh4Cell[cell] = nh4;
        cellVol[cell] = xstep * ystep * zstep;
        ghost[cell] = false;
        if (i == xlo - (xstep/2)) {
        	ghost[cell] = true;
        }
        if (i == xhi + (xstep/2)) {
        	ghost[cell] = true;
        }
        if (j == ylo - (ystep/2)) {
        	ghost[cell] = true;
        }
        if (j == yhi + (ystep/2)) {
        	ghost[cell] = true;
        }
        if (k == zlo - (zstep/2)) {
        	ghost[cell] = true;
        }
        if (k == zhi + (zstep/2)) {
        	ghost[cell] = true;
        }
        cell++;
      }
    }
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

  int cellIn[nall];
  double xHET[numCells];
  double xAOB[numCells];
  double xNOB[numCells];
  double xEPS[numCells];
  double xTot[numCells];

  double R1[numCells];
  double R2[numCells];
  double R3[numCells];
  double R4[numCells];
  double R5[numCells];
  // double R6[numCells] = bHET;
  // double R7[numCells] = bAOB;
  // double R8[numCells] = bNOB;
  // double R9[numCells] = bEPS;
  double Rs[numCells];
  double Ro2[numCells];
  double Rnh4[numCells];
  double Rno2[numCells];
  double Rno3[numCells];
  double cellDo2[numCells];
  double cellDnh4[numCells];
  double cellDno2[numCells];
  double cellDno3[numCells];
  double cellDs[numCells];

  // Figure out which cell each particle is in
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

      int j;
      for (j = 0; j < numCells; j ++) {
        if ((xCell[j] - xstep/2) <= atom->x[i][0] &&
            (xCell[j] + xstep/2) >= atom->x[i][0] &&
            (yCell[j] - ystep/2) <= atom->x[i][1] &&
            (yCell[j] + ystep/2) >= atom->x[i][1] &&
            (zCell[j] - zstep/2) <= atom->x[i][2] &&
            (zCell[j] + zstep/2) >= atom->x[i][2]) {
          cellIn[i] = j;
          xHET[j] += (gHET * rmass[i])/cellVol[j];
          xAOB[j] += (gAOB * rmass[i])/cellVol[j];
          xNOB[j] += (gNOB * rmass[i])/cellVol[j];
          xEPS[j] += (gEPS * rmass[i])/cellVol[j];
          xTot[j] += rmass[i]/cellVol[j];
        }
      }
    }
  }

  for (i = 0; i < numCells; i ++) {
    double diffusionFunction = 1 - ((0.43 * pow(xTot[i], 0.92))/(11.19+0.27*pow(xTot[i], 0.99)));

    cellDo2[i] = diffusionFunction * Do2; 
    cellDnh4[i] = diffusionFunction * Dnh4; 
    cellDno2[i] = diffusionFunction * Dno2; 
    cellDno3[i] = diffusionFunction * Dno3; 
    cellDs[i] = diffusionFunction * Ds; 
  	
    R1[i] = MumHET*(subCell[i]/(KsHET+subCell[i]))*(o2Cell[i]/(Ko2HET+o2Cell[i]));
    R2[i] = MumAOB*(nh4Cell[i]/(Knh4AOB+nh4Cell[i]))*(o2Cell[i]/(Ko2AOB+o2Cell[i]));
    R3[i] = MumNOB*(no2Cell[i]/(Kno2NOB+no2Cell[i]))*(o2Cell[i]/(Ko2NOB+o2Cell[i]));
    R4[i] = etaHET*MumHET*(subCell[i]/(KsHET+subCell[i]))*(no3Cell[i]/(Kno3HET+no3Cell[i]))*(Ko2HET/(Ko2HET+o2Cell[i]));
    R5[i] = etaHET*MumHET*(subCell[i]/(KsHET+subCell[i]))*(no2Cell[i]/(Kno2HET+no2Cell[i]))*(Ko2HET/(Ko2HET+o2Cell[i]));
    Rs[i] = ((-1/YHET)*((R1[i]+R4[i]+R5[i])*xHET[i]))+((1-Y1)*(bHET*xHET[i]+bAOB*xAOB[i]+bNOB*xNOB[i]))+(bEPS*xEPS[i]);
    Ro2[i] = (((1-YHET-YEPS)/YHET)*R1[i]*xHET[i])-(((3.42-YAOB)/YAOB)*R2[i]*xAOB[i])-(((1.15-YNOB)/YNOB)*R3[i]*xNOB[i]);
    Rnh4[i] = -(1/YAOB)*R2[i]*xAOB[i];
    Rno2[i] = ((1/YAOB)*R2[i]*xAOB[i])-((1/YNOB)*R3[i]*xNOB[i])-(((1-YHET-YEPS)/(1.17*YHET))*R5[i]*xHET[i]);
    Rno3[i] = ((1/YNOB)*R3[i]*xNOB[i])-(((1-YHET-YEPS)/(2.86*YHET))*R4[i]*xHET[i]);

    // subCell[i] = subCell[i] - (Rs[i] * update->ntimestep);
    // o2Cell[i] = o2Cell[i] - (Ro2[i] * update->ntimestep);
    // no2Cell[i] = no2Cell[i] - (Rno2[i] * update->ntimestep);
    // no3Cell[i] = no3Cell[i] - (Rno3[i] * update->ntimestep);
    // nh4Cell[i] = nh4Cell[i] - (Rnh4[i] * update->ntimestep);



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

      // int j;
      // for (j = 0; j < numCells; j ++) {
      //   if ((xCell[j] - xstep/2) <= atom->x[i][0] &&
      //       (xCell[j] + xstep/2) >= atom->x[i][0] &&
      //       (yCell[j] - ystep/2) <= atom->x[i][1] &&
      //       (yCell[j] + ystep/2) >= atom->x[i][1] &&
      //       (zCell[j] - zstep/2) <= atom->x[i][2] &&
      //       (zCell[j] + zstep/2) >= atom->x[i][2]) {
      //   	cellIn[i] = j;
      //   	xHET[j] += (gHET * rmass[i])/cellVol[j];
      //   	xAOB[j] += (gAOB * rmass[i])/cellVol[j];
      //   	xNOB[j] += (gNOB * rmass[i])/cellVol[j];
      //   	xEPS[j] += (gEPS * rmass[i])/cellVol[j];
      //    xTot[j] += rmass[i]/cellVol[j];
      //   }
      // }

      // double R1 = MumHET*(sub[i]/(KsHET+sub[i]))*(o2[i]/(Ko2HET+o2[i]));
      // double R2 = MumAOB*(nh4[i]/(Knh4AOB+nh4[i]))*(o2[i]/(Ko2AOB+o2[i]));
      // double R3 = MumNOB*(no2[i]/(Kno2NOB+no2[i]))*(o2[i]/(Ko2NOB+o2[i]));
      // double R4 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no3[i]/(Kno3HET+no3[i]))*(Ko2HET/(Ko2HET+o2[i]));
      // double R5 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no2[i]/(Kno2HET+no2[i]))*(Ko2HET/(Ko2HET+o2[i]));
      // double R6 = bHET;
      // double R7 = bAOB;
      // double R8 = bNOB;
      // double R9 = bEPS;

      double value = update->dt * (gHET*(R1[cellIn[i]]+R4[cellIn[i]]+R5[cellIn[i]]) + gAOB*R2[cellIn[i]] + gNOB*R3[cellIn[i]] - gEPS*bEPS);

      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      double oldMass = rmass[i];
      rmass[i] = rmass[i]*(1 + (value*nevery));
      if (rmass[i] <= 0) {
        rmass[i] = oldMass;
      }

      // fprintf(stdout, "Radius: %e Outer Radius: %e\n", radius[i], outerRadius[i]);
      // fprintf(stdout, "ID: %i Type: %i Outer Mass: %e\n", atom->tag[i], atom->type[i], outerMass[i]);
      

      double value2 = update->dt * (YEPS/YHET)*(R1[cellIn[i]]+R4[cellIn[i]]+R5[cellIn[i]]);
      outerMass[i] = (((4.0*MY_PI/3.0)*((outerRadius[i]*outerRadius[i]*outerRadius[i])-(radius[i]*radius[i]*radius[i])))*EPSdens)+(value2*nevery*rmass[i]);

      outerRadius[i] = pow((3.0/(4.0*MY_PI))*((rmass[i]/density)+(outerMass[i]/EPSdens)),(1.0/3.0));

      radius[i] = pow(((6*rmass[i])/(density*MY_PI)),(1.0/3.0))*0.5;
    }

  }

  modify->addstep_compute(update->ntimestep + nevery);
}
