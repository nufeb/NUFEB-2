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

#include "dump_bio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <dirent.h>

#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "fix.h"
#include "modify.h"
#include "atom.h"
#include "bio.h"
#include "fix_bio_kinetics.h"

#include "compute_bio_diameter.h"
#include "compute_bio_dimension.h"
#include "compute_bio_diversity.h"
#include "compute_bio_height.h"
#include "compute_bio_rough.h"
#include "compute_bio_segregate.h"
#include "compute_bio_ntypes.h"

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
  aveconcFlag = 0;
  catFlag = 0;
  phFlag = 0;
  massFlag = 0;
  massHeader = 0;
  divHeader = 0;
  typeHeader = 0;
  nuHeader = 0;
  gasFlag = 0;
  yieldFlag = 0;

  diaFlag = 0;
  dimFlag = 0;
  divFlag = 0;
  heightFlag = 0;
  roughFlag = 0;
  segFlag = 0;

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
  kinetics = NULL;
  nfix = modify->nfix;
  ncompute = modify->ncompute;

  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required");

  bio = kinetics->bio;

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
    if (strcmp(keywords[i],"concentration") == 0) {
      concFlag = 1;
      if (stat("./Results/S", &st) == -1) {
          mkdir("./Results/S", 0700);
      }
      for (int j = 1; j < nnus + 1; j++) {
        if (bio->nuType[j] == 0 && strcmp(bio->nuName[j], "h") != 0 && strcmp(bio->nuName[j], "h2o") != 0) {
          char *name = bio->nuName[j];
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
        char *name = bio->typeName[j];
        int len = 17;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/DGRAn/");
        strcat(path, name);

        if (stat(path, &st) == -1) {
            mkdir(path, 0700);
        }
      }
    } else if (strcmp(keywords[i],"DGRCat") == 0) {
      catFlag = 1;
      if (stat("./Results/DGRCat", &st) == -1) {
          mkdir("./Results/DGRCat", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        char *name = bio->typeName[j];
        int len = 18;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/DGRCat/");
        strcat(path, name);

        if (stat(path, &st) == -1) {
            mkdir(path, 0700);
        }
      }
    } else if (strcmp(keywords[i],"yield") == 0) {
      yieldFlag = 1;
      if (stat("./Results/yield", &st) == -1) {
          mkdir("./Results/yield", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        char *name = bio->typeName[j];
        int len = 17;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/yield/");
        strcat(path, name);

        if (stat(path, &st) == -1) {
            mkdir(path, 0700);
        }
      }
    } else if (strcmp(keywords[i],"ph") == 0) {
      phFlag = 1;
      if (stat("./Results/pH", &st) == -1) {
          mkdir("./Results/pH", 0700);
      }
    } else if (strcmp(keywords[i],"biomass") == 0) {
      massFlag = 1;
    } else if (strcmp(keywords[i],"ave_concentration") == 0) {
      aveconcFlag = 1;
    } else if (strcmp(keywords[i],"gas") == 0) {
      gasFlag = 1;
      if (stat("./Results/Gas", &st) == -1) {
          mkdir("./Results/Gas", 0700);
      }
      for (int j = 1; j < nnus + 1; j++) {
        if (bio->nuType[j] == 1) {
          char *name = bio->nuName[j];
          int len = 16;
          len += strlen(name);
          char path[len];
          strcpy(path, "./Results/Gas/");
          strcat(path, name);

          if (stat(path, &st) == -1) {
              mkdir(path, 0700);
          }
        }
      }
    } else if (strcmp(keywords[i],"diameter") == 0) {
      diaFlag = 1;
    } else if (strcmp(keywords[i],"dimension") == 0) {
      dimFlag = 1;
    } else if (strcmp(keywords[i],"diversity") == 0) {
      divFlag = 1;
    } else if (strcmp(keywords[i],"ave_height") == 0) {
      heightFlag = 1;
    } else if (strcmp(keywords[i],"roughness") == 0) {
      roughFlag = 1;
    } else if (strcmp(keywords[i],"segregation") == 0) {
      segFlag = 1;
    }  else if (strcmp(keywords[i],"ntypes") == 0) {
      ntypeFlag = 1;
    } else lmp->error->all(FLERR,"Undefined dump_bio keyword");

    i++;
  }

  for (int j = 0; j < ncompute; j++) {
    if (strcmp(modify->compute[j]->style,"diameter") == 0) {
      cdia = static_cast<ComputeNufebDiameter *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"dimension") == 0) {
      cdim = static_cast<ComputeNufebDimension *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"diversity") == 0) {
      cdiv = static_cast<ComputeNufebDiversity *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"ave_height") == 0) {
      cheight = static_cast<ComputeNufebHeight *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"roughness") == 0) {
      crough = static_cast<ComputeNufebRough *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"segregation") == 0) {
      cseg = static_cast<ComputeNufebSegregate *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"ntypes") == 0) {
      ctype = static_cast<ComputeNufebNtypes *>(lmp->modify->compute[j]);
      continue;
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpBio::write()
{
  if (update-> ntimestep == 0) return;

  int nnus = kinetics->bio->nnus;
  int ntypes = atom->ntypes;

  if (concFlag == 1) {
    for (int i = 1; i < nnus+1; i++) {
      if (bio->nuType[i] == 0 && strcmp(bio->nuName[i], "h") != 0 && strcmp(bio->nuName[i], "h2o") != 0) {
        char *name = bio->nuName[i];
        int len = 30;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/S/");
        strcat(path, name);
        strcat(path, "/r*.csv");

        filename = path;
        openfile();
        write_concentration_data(i);
        fclose(fp);
      }
    }
  }

  if (anFlag == 1) {
    for (int i = 1; i < ntypes+1; i++) {
      char *name = bio->typeName[i];
      int len = 50;
      len += strlen(name);
      char path[len];
      strcpy(path, "./Results/DGRAn/");
      strcat(path, name);
      strcat(path, "/r*.csv");

      filename = path;
      openfile();
      write_DGRAn_data(i);
      fclose(fp);
    }
  }

  if (catFlag == 1) {
    for (int i = 1; i < ntypes+1; i++) {
      char *name = bio->typeName[i];
      int len = 50;
      len += strlen(name);
      char path[len];
      strcpy(path, "./Results/DGRCat/");
      strcat(path, name);
      strcat(path, "/r*.csv");

      filename = path;
      openfile();
      write_DGRCat_data(i);
      fclose(fp);
    }
  }

  if (yieldFlag == 1) {
    for (int i = 1; i < ntypes+1; i++) {
      char *name = bio->typeName[i];
      int len = 50;
      len += strlen(name);
      char path[len];
      strcpy(path, "./Results/yield/");
      strcat(path, name);
      strcat(path, "/r*.csv");

      filename = path;
      openfile();
      write_yield_data(i);
      fclose(fp);
    }
  }

  if (phFlag == 1) {
    int len = 30;
    char path[len];
    strcpy(path, "./Results/pH/r*.csv");

    filename = path;
    openfile();
    write_pH_data();
    fclose(fp);
  }

  if (massFlag == 1) {
    int len = 35;
    char path[len];
    strcpy(path, "./Results/biomass.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_biomass_data();
    fclose(fp);
  }

  if (diaFlag == 1) {
    int len = 36;
    char path[len];
    strcpy(path, "./Results/diameter.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_diameter_data();
    fclose(fp);
  }

  if (dimFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/dimension.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_dimension_data();
    fclose(fp);
  }

  if (divFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/diversity.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_diversity_data();
    fclose(fp);
  }

  if (heightFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ave_height.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_height_data();
    fclose(fp);
  }

  if (roughFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/roughness.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_rough_data();
    fclose(fp);
  }

  if (segFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/segregation.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_segregate_data();
    fclose(fp);
  }

  if (ntypeFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ntypes.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_ntype_data();
    fclose(fp);
  }

  if (aveconcFlag == 1) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ave_concentration.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_aveconcentration_data();
    fclose(fp);
  }

  if (gasFlag == 1) {
    for (int i = 1; i < nnus+1; i++) {
      if (bio->nuType[i] == 1) {
        char *name = bio->nuName[i];
        int len = 30;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/Gas/");
        strcat(path, name);
        strcat(path, "/r*.csv");

        filename = path;
        openfile();
        write_gas_data(i);
        fclose(fp);
      }
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
  if (nprocs > 1)
  {
    sprintf(filecurrent,"%s_%d_" BIGINT_FORMAT "%s",
            filestar,me,update->ntimestep,ptr+1);
  }
  else
  {
    sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
            filestar,update->ntimestep,ptr+1);
  }
  *ptr = '*';
  fp = fopen(filecurrent,"w");
  delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpBio::pack(tagint *ids)
{
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_concentration_data(int nuID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%e\n",i, x, y, z, kinetics->nuS[nuID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_DGRCat_data(int typeID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%e\n",i, x, y, z, kinetics->DRGCat[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_yield_data(int typeID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%e\n",i, x, y, z, kinetics->gYield[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_DGRAn_data(int typeID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%e\n",i, x, y, z, kinetics->DRGAn[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_pH_data()
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%e\n",i, x, y, z, -log10(kinetics->Sh[i]));
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_aveconcentration_data()
{
  if (!nuHeader) {
    for(int i = 1; i < kinetics->nnus+1; i++){
      fprintf(fp, "%s,\t", kinetics->bio->nuName[i]);
    }
    nuHeader = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%i,\t", update->ntimestep);

  for(int i = 1; i < kinetics->nnus+1; i++){
    double s = 0;
    for(int j = 0; j < kinetics->bgrids; j++){
      s += kinetics->nuS[i][j];
    }
    s = s/kinetics->bgrids;
    fprintf(fp, "%e,\t", s);
  }
  fprintf(fp, "\n");
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_biomass_data()
{
  if (!massHeader) {
    for(int i = 1; i < atom->ntypes+1; i++){
      fprintf(fp, "%s,\t", kinetics->bio->typeName[i]);
    }
    massHeader = 1;
    fprintf(fp, "\n");
  }

  int local = atom->nlocal;
  int ghost = atom->nghost;
  int all = local + ghost;
  double *mass = new double[atom->ntypes+1]();

  for(int i = 0; i < all; i++){
    int type = atom->type[i];
    mass[type] += atom->rmass[i];
  }

  fprintf(fp, "%i,\t", update->ntimestep);

  for(int i = 1; i < atom->ntypes+1; i++){
    fprintf(fp, "%e,\t", mass[i]);
  }
  fprintf(fp, "\n");

  if (update->ntimestep == update->nsteps) {
    for(int i = 1; i < atom->ntypes+1; i++){
      fprintf(fp, "%s,\t", kinetics->bio->typeName[i]);
    }
    massHeader = 1;
    fprintf(fp, "\n");
  }

  delete[] mass;
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_diameter_data()
{
  cdia->compute_scalar();
  fprintf(fp, "%i,\t %e, \n", update->ntimestep, cdia->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_dimension_data()
{
  cdim->compute_scalar();
  fprintf(fp, "%i,\t %e, \n", update->ntimestep, cdim->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_diversity_data()
{
  cdiv->compute_scalar();
  fprintf(fp, "%i,\t %e, \n", update->ntimestep, cdiv->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_height_data()
{
  cheight->compute_scalar();
  fprintf(fp, "%i,\t %e, \n", update->ntimestep, cheight->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_rough_data()
{
  crough->compute_scalar();
  fprintf(fp, "%i,\t %e, \n", update->ntimestep, crough->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_segregate_data()
{
  cseg->compute_scalar();
  fprintf(fp, "%i,\t %e, \n", update->ntimestep, cseg->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_ntype_data()
{
  if (!typeHeader) {
    for(int i = 1; i < atom->ntypes+1; i++){
      fprintf(fp, "%s,\t", kinetics->bio->typeName[i]);
    }
    typeHeader = 1;
    fprintf(fp, "\n");
  }

  int local = atom->nlocal;
  int ghost = atom->nghost;
  int all = local + ghost;
  int *ntypes1 = new int[atom->ntypes+1]();

  for(int i = 0; i < all; i++){
    int type = atom->type[i];
    ntypes1[type] ++;
  }

  fprintf(fp, "%i,\t", update->ntimestep);

  for(int i = 1; i < atom->ntypes+1; i++){
    fprintf(fp, "%i,\t", ntypes1[i]);
  }
  fprintf(fp, "\n");

  delete[] ntypes1;
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_gas_data(int nuID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    fprintf(fp, "%i,\t%f,\t%f,\t%f,\t%e\n",i, x, y, z, kinetics->qGas[nuID][i]);
  }
}


bigint DumpBio::memory_usage() {
  return 0;
}

