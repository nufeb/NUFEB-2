/*
 * fix_metabolism.h
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

#include "fix_kinetics_thermo.h"

#include <math.h>
#include <string.h>
#include <cstdio>
#include <string>
#include <sstream>

#include "atom.h"
#include "bio.h"
#include "domain.h"
#include "error.h"
#include "fix_kinetics.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "memory.h"
#include "modify.h"
#include "pointers.h"
#include "variable.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;


/* ---------------------------------------------------------------------- */

FixKineticsThermo::FixKineticsThermo(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Not enough arguments in fix kinetics/thermo command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix kinetics/thermo command");
}

/* ---------------------------------------------------------------------- */

FixKineticsThermo::~FixKineticsThermo()
{
  memory->destroy(dG0);
  memory->destroy(DRGCat);
  memory->destroy(DRGAn);
  memory->destroy(khV);
}

/* ---------------------------------------------------------------------- */

int FixKineticsThermo::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::init()
{
  // register fix kinetics with this class
  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  ntypes = atom->ntypes;

  bio = kinetics->bio;
  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;
  ngrids = nx * ny * nz;
  temp = kinetics->temp;
  rth = kinetics->rth;
  metCoeff = kinetics->metCoeff;

  nnus = bio->nnus;
  catCoeff = bio->catCoeff;
  anabCoeff = bio->anabCoeff;
  nuGCoeff = bio->nuGCoeff;
  typeGCoeff = bio->typeGCoeff;
  yield = bio->yield;
  diss = bio->dissipation;
  ngflag = bio->ngflag;
  tgflag = bio->tgflag;

  DRGCat = memory->create(DRGCat,ntypes+1,ngrids,"kinetics/thermo:DRGCat");
  DRGAn = memory->create(DRGAn,ntypes+1,ngrids,"kinetics/thermo:DRGAn");
  khV = memory->create(khV,nnus+1,"kinetics/thermo:khV");

  init_dG0();
  init_KhV();
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_dG0()
{
  dG0 = memory->create(dG0,ntypes+1,2,"kinetics/thermo:dG0");

  for (int i = 1; i <= ntypes; i++) {
    dG0[i][0] = 0;
    dG0[i][1] = 0;

    for (int j = 1; j <= nnus; j++) {
      int flag = ngflag[j];
      dG0[i][0] += catCoeff[i][j] * nuGCoeff[j][flag];
      dG0[i][1] += anabCoeff[i][j] * nuGCoeff[j][flag];
    }

    int flag = tgflag[i];
    dG0[i][1] += typeGCoeff[i][flag];
//    printf("dG0[%i][0] = %e \n", i, dG0[i][0]);
//    printf("dG0[%i][1] = %e \n", i, dG0[i][1]);
  }
}

/* ----------------------------------------------------------------------*/

void FixKineticsThermo::init_KhV()
{
  for (int i = 1; i < nnus + 1; i++) {
    khV[i] = 0;
  }

  for (int i = 1; i < nnus + 1; i++) {
    char *nName;     // nutrient name
    nName = bio->nuName[i];

    if (bio->nuType[i] == 1) {
      char *lName = new char[strlen(nName)];     // corresponding liquid
      strncpy(lName, nName+1, strlen(nName));

      for (int j = 1; j < nnus + 1; j++) {
        char *nName2;     // nutrient name
        nName2 = bio->nuName[j];
        if (strcmp(lName, nName2) == 0) {
          if (nuGCoeff[j][0] > 10000) {
            khV[j] = exp((nuGCoeff[j][1] - nuGCoeff[i][1]) / (-rth * temp));
          } else {
            khV[j] = exp((nuGCoeff[j][0] - nuGCoeff[i][1]) / (-rth * temp));
          }
          break;
        }
      }
      memory->sfree(lName);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixKineticsThermo::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  thermo();
  output_data();
}

/* ----------------------------------------------------------------------
  thermodynamics
------------------------------------------------------------------------- */
void FixKineticsThermo::thermo()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int *type = atom->type;

  iyield = kinetics->iyield;
  nuS = kinetics->nuS;

  for (int i = 0; i < ngrids; i++) {
    for (int j = 1; j <= ntypes; j++) {
      //Gibbs free energy of the reaction
      DRGCat[j][i] = dG0[j][0];  //catabolic energy values
      DRGAn[j][i] = dG0[j][1];  //anabolic energy values

      for (int k = 1; k <= nnus; k++) {
        if (strcmp(bio->nuName[k], "h2o") != 0){
          if (nuGCoeff[k][1] < 1e4) {
            double x = temp * rth * log(nuS[k][i]);
            DRGCat[j][i] += catCoeff[j][k] * x;
            DRGAn[j][i] += anabCoeff[j][k] * x;
          }
        }
      }
//      printf("dGrxn[%i][%i][0] = %e \n", i,j, DRGAn[j][i]);
//      printf("dGrxn[%i][%i][1] = %e \n", i,j, DRGCat[j][i]);

      //use catabolic and anabolic energy values to derive catabolic reaction equation
      if (DRGCat[j][i] < 0)
        iyield[j][i] = -(DRGAn[j][i] + diss[j]) / DRGCat[j][i];
      else
        iyield[j][i] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
  output energy to data file
------------------------------------------------------------------------- */

void FixKineticsThermo::output_data(){
  std::ostringstream stm;
  stm << update->ntimestep;
  string str = "./DGRCat/DGRCat.csv."+ stm.str();
  FILE *pFile = fopen (str.c_str(), "a");
  fprintf(pFile, ",x,y,z,scalar,1,1,1,0.5\n");
  double average = 0.0;

  double xlo,xhi,ylo,yhi,zlo,zhi;

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

  double stepx = (xhi - xlo) / nx;
  double stepy = (yhi - ylo) / ny;
  double stepz = (zhi - zlo) / nz;

  for(int i = 0; i < ngrids; i++){
    int zpos = i/(nx * ny) + 1;
    int ypos = (i - (zpos - 1) * (nx * ny)) / nx + 1;
    int xpos = i - (zpos - 1) * (nx * ny) - (ypos - 1) * nx + 1;

    double x = xpos * stepx - stepx/2;
    double y = ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    average += DRGCat[2][i];

    fprintf(pFile, "%i,\t%f,\t%f,\t%f,\t%f\n",i, x, y, z, DRGCat[2][i]);
  }
  fclose(pFile);

  average = average / ngrids;
  string str2 = "./DGRCat/Average.csv";
  pFile = fopen (str2.c_str(), "a");
  fprintf(pFile, "%e\n", average);

}
