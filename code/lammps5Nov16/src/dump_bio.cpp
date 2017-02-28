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
#include "dump_bio.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "fix.h"

using namespace LAMMPS_NS;

// customize by adding keyword

#define NOUTPUTS 2

/* ---------------------------------------------------------------------- */

DumpBio::DumpBio(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg == 4) error->all(FLERR,"No dump bio arguments specified");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal dump bio command");

  nfix = 0;
  id_fix = NULL;
  fix = NULL;
  selected_outputs = NULL;
  fp = NULL;

  // customize for new sections
  selected_outputs = (char **) memory->srealloc(selected_outputs, 2*sizeof(char *), "selection");

  for (int i = 0; i < 2; i++) {
    int n = strlen(arg[i+4]) + 3;
    selected_outputs[i] = new char[n];
    strcpy(selected_outputs[i], arg[i+4]);
    strcat(selected_outputs[i],".*");
  }

//  for (int i = 0; i < 2; i++) {
//    printf("%s \n", selected_outputs[i]);
//  }
}

/* ---------------------------------------------------------------------- */

DumpBio::~DumpBio()
{
  for (int i = 0; i < 2; i++) {
    delete [] selected_outputs[i];
  }

  memory->sfree(selected_outputs);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write()
{
  for (int i = 0; i < 2; i++) {
    filename = selected_outputs[i];
   // printf("%s \n", filename);
    openfile();
    //fprintf(fp,"ITEM: ATOMS\n");
    fclose(fp);
  }
  // if file per timestep, close file if I am filewriter
}

/* ---------------------------------------------------------------------- */
void DumpBio::openfile()
{

  // if one file per timestep, replace '*' with current timestep
  char *filecurrent = filename;
  printf("%s \n", filename);
  char *filestar = filecurrent;
  filecurrent = new char[strlen(filestar) + 16];
  char *ptr = strchr(filestar,'*');
  *ptr = '\0';

  sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
          filestar,update->ntimestep,ptr+1);
  *ptr = '*';

  fp = fopen(filecurrent,"w");
  printf("here! \n");
  delete [] filecurrent;
}


/* ---------------------------------------------------------------------- */

void DumpBio::write_header(bigint ndump)
{
  //fprintf(fp,"ITEM: ATOMS\n");
}

/* ---------------------------------------------------------------------- */
void DumpBio::init_style()
{

}

/* ---------------------------------------------------------------------- */

void DumpBio::write_data(int n, double *mybuf)
{
  //fprintf(fp,"ITEM: ATOMS\n");
}

/* ---------------------------------------------------------------------- */

void DumpBio::pack(tagint *ids)
{
}
