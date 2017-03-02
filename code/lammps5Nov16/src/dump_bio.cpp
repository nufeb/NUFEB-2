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

/* ----------------------------------------------------------------------
   Contributing author:  Bowen LI
------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "dump_bio.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "fix.h"
#include "fix_kinetics.h"
#include "modify.h"
#include "bio.h"
#include "atom.h"

struct stat st = {0};

using namespace LAMMPS_NS;

// customize by adding keyword

/* ---------------------------------------------------------------------- */

DumpBio::DumpBio(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg <= 4) error->all(FLERR,"No dump bio arguments specified");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal dump bio command");

  nkeywords = narg - 4;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;
  keywords = NULL;
  fp = NULL;

  anFlag = 0;
  concFlag = 0;
  catFlag = 0;
  phFlag = 0;

  // customize for new sections
  keywords = (char **) memory->srealloc(keywords, nkeywords*sizeof(char *), "keywords");

  for (int i = 0; i < nkeywords; i++) {
    int n = strlen(arg[i+4]) + 2;
    keywords[i] = new char[n];
    strcpy(keywords[i], arg[i+4]);
  }
}

/* ---------------------------------------------------------------------- */

DumpBio::~DumpBio()
{
  for (int i = 0; i < nkeywords; i++) {
    delete [] keywords[i];
  }

  memory->sfree(keywords);
}


/* ---------------------------------------------------------------------- */
void DumpBio::init_style()
{
  // register fix kinetics with this class
  nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  ptr = kinetics->bio;

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
  }
  else {
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

  int i = 0;
  int nnus = kinetics->bio->nnus;
  int ntypes = atom->ntypes;

  //create directory
  if (stat("./Results", &st) == -1) {
      mkdir("./Results", 0700);
  }

  while (i < nkeywords) {
    if (strcmp(keywords[i],"Conc") == 0) {
      concFlag = 1;
      if (stat("./Results/S", &st) == -1) {
          mkdir("./Results/S", 0700);
      }
      for (int j = 1; j < nnus + 1; j++) {
        if (ptr->nuType[j] == 0 && strcmp(ptr->nuName[j], "h") != 0 && strcmp(ptr->nuName[j], "h2o") != 0) {
          char *name = ptr->nuName[j];
          int len = 13;
          len += strlen(name);
          char path[len];
          strcpy(path, "./Results/S/");
          strcat(path, name);

          if (stat(path, &st) == -1) {
              mkdir(path, 0700);
          }
        }
      }
    } else if (strcmp(keywords[i],"DGRAn") == 0) {
      anFlag = 1;
      if (stat("./Results/DGRAn", &st) == -1) {
          mkdir("./Results/DGRAn", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        if (ptr->nuType[j] == 0) {
          char *name = ptr->typeName[j];
          int len = 17;
          len += strlen(name);
          char path[len];
          strcpy(path, "./Results/DGRAn/");
          strcat(path, name);

          if (stat(path, &st) == -1) {
              mkdir(path, 0700);
          }
        }
      }
    } else if (strcmp(keywords[i],"DGRCat") == 0) {
      catFlag = 1;
      if (stat("./Results/DGRCat", &st) == -1) {
          mkdir("./Results/DGRCat", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        if (ptr->nuType[j] == 0) {
          char *name = ptr->typeName[j];
          int len = 18;
          len += strlen(name);
          char path[len];
          strcpy(path, "./Results/DGRCat/");
          strcat(path, name);

          if (stat(path, &st) == -1) {
              mkdir(path, 0700);
          }
        }
      }
    }

    i++;
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write()
{
  int i = 0;
  int nnus = kinetics->bio->nnus;
  int ntypes = atom->ntypes;

  if (concFlag == 1) {
    for (int j = 1; j < nnus+1; j++) {
      if (ptr->nuType[j] == 0 && strcmp(ptr->nuName[j], "h") != 0 && strcmp(ptr->nuName[j], "h2o") != 0) {
        char *name = ptr->nuName[j];
        int len = 30;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/S/");
        strcat(path, name);
        strcat(path, "/*.csv");

        filename = path;
        openfile();
        write_diffsuion_data(j);
        fclose(fp);
      }
    }
  }

  if (anFlag == 1) {
    for (int j = 1; j < ntypes+1; j++) {
      char *name = ptr->typeName[j];
      int len = 50;
      len += strlen(name);
      char path[len];
      strcpy(path, "./Results/DGRAn/");
      strcat(path, name);
      strcat(path, "/*.csv");

      filename = path;
      openfile();
      write_DGRAn_data(j);
      fclose(fp);
    }
  }

  if (catFlag == 1) {
    for (int j = 1; j < ntypes+1; j++) {
      char *name = ptr->typeName[j];
      int len = 50;
      len += strlen(name);
      char path[len];
      strcpy(path, "./Results/DGRCat/");
      strcat(path, name);
      strcat(path, "/*.csv");

      filename = path;
      openfile();
      write_DGRCat_data(j);
      fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_header(bigint n)
{
}

void DumpBio::write_data(int n, double *mybuf)
{
}
/* ---------------------------------------------------------------------- */
void DumpBio::openfile()
{
  //replace '*' with current timestep
  char *filecurrent = filename;
  //printf("%s \n", filename);
  char *filestar = filecurrent;
  filecurrent = new char[strlen(filestar) + 16];
  char *ptr = strchr(filestar,'*');
  *ptr = '\0';

  sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
          filestar,update->ntimestep,ptr+1);
  *ptr = '*';
  fp = fopen(filecurrent,"w");
  printf("%s, \n", filecurrent);
  delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpBio::pack(tagint *ids)
{
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_diffsuion_data(int nuID)
{
  fprintf(fp, "x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(nx * ny) + 1;
    int ypos = (i - (zpos - 1) * (nx * ny)) / nx + 1;
    int xpos = i - (zpos - 1) * (nx * ny) - (ypos - 1) * nx + 1;

    double x = xpos * stepx - stepx/2;
    double y = ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%f\n",i, x, y, z, kinetics->nuS[nuID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_DGRCat_data(int typeID)
{
  fprintf(fp, "x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(nx * ny) + 1;
    int ypos = (i - (zpos - 1) * (nx * ny)) / nx + 1;
    int xpos = i - (zpos - 1) * (nx * ny) - (ypos - 1) * nx + 1;

    double x = xpos * stepx - stepx/2;
    double y = ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%f\n",i, x, y, z, kinetics->DRGCat[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_DGRAn_data(int typeID)
{
  fprintf(fp, "x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(nx * ny) + 1;
    int ypos = (i - (zpos - 1) * (nx * ny)) / nx + 1;
    int xpos = i - (zpos - 1) * (nx * ny) - (ypos - 1) * nx + 1;

    double x = xpos * stepx - stepx/2;
    double y = ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%f\n",i, x, y, z, kinetics->DRGAn[typeID][i]);
  }
}

bigint DumpBio::memory_usage() {
  return 0;
}
