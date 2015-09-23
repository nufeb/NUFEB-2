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
  if (narg < 28) error->all(FLERR,"Illegal fix growth command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix growth command");

  var = new char*[21];
  ivar = new int[21];

  int i;
  for (i = 0; i < 21; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  nx = atoi(arg[25]);
  ny = atoi(arg[26]);
  nz = atoi(arg[27]);

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

  xloBound = false;
  xhiBound = false;
  yloBound = false;
  yhiBound = false;
  zloBound = false;
  zhiBound = false;

  if (narg > 28) {
  	for (i = 29; i < narg; i++) {
  		if (strcmp(arg[i],"xlo") == 0) {
  			xloBound = true;
  		}
  		else if (strcmp(arg[i],"xhi") == 0) {
  			xhiBound = true;
  		}
  		else if (strcmp(arg[i],"ylo") == 0) {
  			yloBound = true;
  		}
  		else if (strcmp(arg[i],"yhi") == 0) {
  			yhiBound = true;
  		}
  		else if (strcmp(arg[i],"zlo") == 0) {
  			zloBound = true;
  		}
  		else if (strcmp(arg[i],"zhi") == 0) {
  			zhiBound = true;
  		}
  		else {
  			error->all(FLERR,"Boundary must be xlo, xhi, ylo, yhi, zlo, or zhi.");
  		}
  	}
  }

}

/* ---------------------------------------------------------------------- */

FixDiffNuGrowth::~FixDiffNuGrowth()
{
  int i;
  for (i = 0; i < 21; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
  delete [] xCell;
  delete [] yCell;
  delete [] zCell;
  delete [] cellVol;
  delete [] boundary;
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
  for (n = 0; n < 21; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR,"Variable name for fix nugrowth does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR,"Variable for fix nugrowth is invalid style");
  }

  numCells = nx*ny*nz;

  xCell = new double[numCells];
  yCell = new double[numCells];
  zCell = new double[numCells];
  cellVol = new double[numCells];
  boundary = new bool[numCells];
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
  for (i = xlo + (xstep/2); i < xhi; i += xstep) {
    for (j = ylo + (ystep/2); j < yhi; j += ystep) {
      for (k = zlo + (zstep/2); k < zhi; k += zstep) {
        xCell[cell] = i;
        yCell[cell] = j;
        zCell[cell] = k;
        cellVol[cell] = xstep * ystep * zstep;
        boundary[cell] = false;
        if (i == xlo + (xstep/2) && xloBound) {
        	boundary[cell] = true;
        }
        if (i == xhi - (xstep/2) && xhiBound) {
        	boundary[cell] = true;
        }
        if (j == ylo + (ystep/2) && yloBound) {
        	boundary[cell] = true;
        }
        if (j == yhi - (ystep/2) && yhiBound) {
        	boundary[cell] = true;
        }
        if (k == zlo + (zstep/2) && zloBound) {
        	boundary[cell] = true;
        }
        if (k == zhi - (zstep/2) && zhiBound) {
        	boundary[cell] = true;
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
  // double bHET = input->variable->compute_equal(ivar[12]);
  // double bAOB = input->variable->compute_equal(ivar[13]);
  // double bNOB = input->variable->compute_equal(ivar[14]);
  double bEPS = input->variable->compute_equal(ivar[12]);
  double YEPS = input->variable->compute_equal(ivar[13]);
  double YHET = input->variable->compute_equal(ivar[14]);
  double EPSdens = input->variable->compute_equal(ivar[15]);
  double Do2 = input->variable->compute_equal(ivar[16]);
  double Dnh4 = input->variable->compute_equal(ivar[17]);
  double Dno2 = input->variable->compute_equal(ivar[18]);
  double Dno3 = input->variable->compute_equal(ivar[19]);
  double Ds = input->variable->compute_equal(ivar[20]);

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
  double cHET[numCells];
  double cAOB[numCells];
  double cNOB[numCells];
  double cEPS[numCells];

  double R1[numCells];
  double R2[numCells];
  double R3[numCells];
  double R4[numCells];
  double R5[numCells];
  double R6[numCells];
  double R7[numCells];
  double R8[numCells];
  double R9[numCells];
  double Ro2[numCells];
  double Rnh4[numCells];
  double Rno2[numCells];
  double Rno3[numCells];
  double Rs[numCells];

  // Figure out which cell each particle is in

  for (i = 0; i < numCells; i ++) {
  	// Calculate R's at the cell level using the nutrients
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

      int j;
      for (j = 0; j < numCells; j ++) {
        if ((xCell[j] - xstep/2) <= atom->x[i][0] &&
            (xCell[j] + xstep/2) >= atom->x[i][0] &&
            (yCell[j] - ystep/2) <= atom->x[i][1] &&
            (yCell[j] + ystep/2) >= atom->x[i][1] &&
            (zCell[j] - zstep/2) <= atom->x[i][2] &&
            (zCell[j] + zstep/2) >= atom->x[i][2]) {
        	cellIn[i] = j;
        	cHET[j] += (gHET * rmass[i])/cellVol[j];
        	cAOB[j] += (gAOB * rmass[i])/cellVol[j];
        	cNOB[j] += (gNOB * rmass[i])/cellVol[j];
        	cEPS[j] += (gEPS * rmass[i])/cellVol[j];
        }
      }

      double R1 = MumHET*(sub[i]/(KsHET+sub[i]))*(o2[i]/(Ko2HET+o2[i]));
      double R2 = MumAOB*(nh4[i]/(Knh4AOB+nh4[i]))*(o2[i]/(Ko2AOB+o2[i]));
      double R3 = MumNOB*(no2[i]/(Kno2NOB+no2[i]))*(o2[i]/(Ko2NOB+o2[i]));
      double R4 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no3[i]/(Kno3HET+no3[i]))*(Ko2HET/(Ko2HET+o2[i]));
      double R5 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no2[i]/(Kno2HET+no2[i]))*(Ko2HET/(Ko2HET+o2[i]));
      // double R6 = bHET;
      // double R7 = bAOB;
      // double R8 = bNOB;
      double R9 = bEPS;

      double value = update->dt * (gHET*(R1+R4+R5) + gAOB*R2 + gNOB*R3 - gEPS*R9);

      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      double oldMass = rmass[i];
      rmass[i] = rmass[i]*(1 + (value*nevery));
      if (rmass[i] <= 0) {
        rmass[i] = oldMass;
      }

      // fprintf(stdout, "Radius: %e Outer Radius: %e\n", radius[i], outerRadius[i]);
      // fprintf(stdout, "ID: %i Type: %i Outer Mass: %e\n", atom->tag[i], atom->type[i], outerMass[i]);
      

      double value2 = update->dt * (YEPS/YHET)*(R1+R4+R5);
      outerMass[i] = (((4.0*MY_PI/3.0)*((outerRadius[i]*outerRadius[i]*outerRadius[i])-(radius[i]*radius[i]*radius[i])))*EPSdens)+(value2*nevery*rmass[i]);

      outerRadius[i] = pow((3.0/(4.0*MY_PI))*((rmass[i]/density)+(outerMass[i]/EPSdens)),(1.0/3.0));

      radius[i] = pow(((6*rmass[i])/(density*MY_PI)),(1.0/3.0))*0.5;
    }

  }

  modify->addstep_compute(update->ntimestep + nevery);
}
