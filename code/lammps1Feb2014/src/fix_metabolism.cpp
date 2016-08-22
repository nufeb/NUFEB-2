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

#include "fix_metabolism.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
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
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define CHUNK 1000
#define MAXLINE 256

#define NSECTIONS 4

/* ---------------------------------------------------------------------- */

FixMetabolism::FixMetabolism(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];

  readfile(arg[3]);

}

/* ---------------------------------------------------------------------- */

FixMetabolism::~FixMetabolism()
{
  memory->destroy(cellConc);
  memory->destroy(bcConc);
  memory->destroy(diffCoeff);

  for (int i = 0; i < atom->ntypes+1; i++) {
    delete [] catCoeff[i];
    delete [] anabCoeff[i];
  }

  memory->sfree(catCoeff);
  memory->sfree(anabCoeff);
  //memory->sfree(metaCoeff);

  for (int i = 0; i < nsubs+1; i++) {
    delete [] subName[i];
  }

  memory->sfree(subName);
}

/* ---------------------------------------------------------------------- */

int FixMetabolism::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMetabolism::init()
{

}

/* ---------------------------------------------------------------------- */

void FixMetabolism::pre_force(int vflag)
{

}

/* ---------------------------------------------------------------------- */

void FixMetabolism::change_dia()
{

}

/* ----------------------------------------------------------------------
   read target coordinates from file, store with appropriate atom
------------------------------------------------------------------------- */

void FixMetabolism::readfile(char *file)
{
  int firstpass = 1;
  int subflag;

  subflag = 0;
  while (1) {
    // open file on proc 0

    if (me == 0) {
      if (firstpass && screen) fprintf(screen,"Reading subtrate file ...\n");
      open(file);
    } else fp = NULL;

    // read header info

    header();

    // allocate arrays

    if (firstpass) {
      allocate_arrays();
    }

    while (strlen(keyword)) {

      if (strcmp(keyword,"Substrates") == 0) {
        subflag = 1;
        if (firstpass) substrate();
        else skip_lines(nsubs);
      } else if (strcmp(keyword,"Diffusion") == 0) {
        if (subflag == 0) error->all(FLERR,"Must read Substrates before Lines");
        if (firstpass) diffusion();
        else skip_lines(nsubs);
      } else if (strcmp(keyword,"Cat") == 0) {
        if (subflag == 0) error->all(FLERR,"Must read Substrates before Lines");
        if (firstpass) cat();
        else skip_lines(atom->ntypes);
      } else if (strcmp(keyword,"Anab") == 0) {
        if (subflag == 0) error->all(FLERR,"Must read Substrates before Lines");
        if (firstpass) anab();
        else skip_lines(atom->ntypes);
    }  else {
      char str[128];
      sprintf(str,"Unknown identifier in data file: %s",keyword);
      error->all(FLERR,str);
    }

    parse_keyword(0);
    }

    // error if natoms > 0 yet no atoms were read

    if (nsubs > 0 && subflag == 0)
      error->all(FLERR,"No Substrates in data file");

    // close file

    if (me == 0) {
      if (compressed) pclose(fp);
      else fclose(fp);
    }

    // done if this was 2nd pass

    if (!firstpass) break;
  }
}


/*------------------------------------------------------------------------- */

void FixMetabolism::header()
{
  int n;
  char *ptr;

  // customize for new sections

  const char *section_keywords[NSECTIONS] =
    {"Substrates", "Diffusion", "Cat", "Anab"};

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
  }
  while (1) {

    // read a line and bcast length if flag is set

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) n = 0;
      else n = strlen(line) + 1;
    }


    MPI_Bcast(&n,1,MPI_INT,0,world);

    // if n = 0 then end-of-file so return with blank line

    if (n == 0) {
      line[0] = '\0';
      return;
    }
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // trim anything from '#' onward
    // if line is blank, continue

    if (ptr = strchr(line,'#')) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keyword and set corresponding variable

    if (strstr(line,"substrates")) {
      sscanf(line,"%i",&nsubs);
    } else break;
  }

  // check that exiting string is a valid section keyword

  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword,section_keywords[n]) == 0) break;
  if (n == NSECTIONS) {
    char str[128];
    sprintf(str,"Unknown identifier in data file: %s",keyword);
    error->all(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   proc 0 opens substrate data file
   test if gzipped
------------------------------------------------------------------------- */

void FixMetabolism::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix,".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file,"r");
  else {
#ifdef LAMMPS_GZIP
    char gunzip[128];
    sprintf(gunzip,"gzip -c -d %s",file);

#ifdef _WIN32
    fp = _popen(gunzip,"rb");
#else
    fp = popen(gunzip,"r");
#endif

#else
    error->one(FLERR,"Cannot open gzipped file");
#endif
  }

  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open file %s",file);
    error->one(FLERR,str);
  }
}

void FixMetabolism::allocate_arrays()
{
  subName = (char **) memory->srealloc(atom->typeName,(nsubs+1)*sizeof(char *),
                                     "atom:typeName");
  cellConc = new double[nsubs+1];
  bcConc = new double[nsubs+1];
  diffCoeff = new double[nsubs+1];

  for (int i = 0; i < atom->ntypes+1; i++) {
    catCoeff[i] = new int[nsubs+1];
    anabCoeff[i] = new int[nsubs+1];
  }

}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
------------------------------------------------------------------------- */

void FixMetabolism::parse_keyword(int first)
{
  int eof = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == NULL) eof = 1;
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof,1,MPI_INT,0,world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   could be skipping Natoms lines, so use bigints
------------------------------------------------------------------------- */

void FixMetabolism::skip_lines(bigint n)
{
  if (me) return;
  char *eof;
  for (bigint i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == NULL) error->one(FLERR,"Unexpected end of data file");
}

void FixMetabolism::diffusion()
{
  int i,m;
  char *next;
  char *buf = new char[atom->ntypes*MAXLINE];

  int eof = comm->read_lines_from_file(fp,nsubs,MAXLINE,buf);
  if (eof) error->all(FLERR,"Unexpected end of data file");

  char *original = buf;
  for (i = 0; i < nsubs; i++) {
    next = strchr(buf,'\n');
    *next = '\0';
    atom->set_ks(buf);
    buf = next + 1;
  }
  delete [] original;
}
//
//void FixMetabolism::set_diffusion(const char *str)
//{
//  if (diffCoeff == NULL) error->all(FLERR,"Cannot set ks value for this atom style");
//
//  char* typeName;
//  double diffCoeff_one;
//  int len = strlen(str);
//  typeName = new char[len];
//
//  int n = sscanf(str,"%s %lg",typeName,&diffCoeff_one);
//  if (n != 2) error->all(FLERR,"Invalid growth line in data file");
//
//  int itype = find_typeID(typeName);
//  if (itype < 1 || itype > ntypes)
//    error->all(FLERR,"Invalid type for diffCoeff set");
//
//  diffCoeff[itype] = diffCoeff_one;
//  //mass_setflag[itype] = 1;
//
//  if (ks[itype] <= 0.0) error->all(FLERR,"Invalid growth value");
//}
//
//int FixMetabolism::find_subID(char *name) {
//
//  for (int i = 0; i < atom->ntypes+1; i++)
//    if (typeName[i] && strcmp(typeName[i],name) == 0) {
//      return i;
//    }
//
//  return -1;
//}
