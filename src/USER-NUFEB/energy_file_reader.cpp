/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#include "energy_file_reader.h"

#include <string.h>
#include "error.h"
#include "atom.h"
#include "memory.h"
#include "comm.h"
#include "grid.h"

using namespace LAMMPS_NS;

#define MAXLINE 256
#define CHUNK 1024
#define DELTA 4            // must be 2 or larger
#define NSECTIONS 12       // change when add to header::section_keywords

EnergyFileReader::EnergyFileReader(LAMMPS *lmp, char *filename) :
        Pointers(lmp),
        filename(filename) {
  MPI_Comm_rank(world,&me);
  ngroups = 0;

  line = new char[MAXLINE];
  keyword = new char[MAXLINE];
  style = new char[MAXLINE];
  buffer = new char[CHUNK*MAXLINE];
  fp = nullptr;
  narg = maxarg = 0;
  arg = nullptr;

  uptake = 0.0;
  max_yield = 1.0;
  decay = 0.0;
  maintain = 0.0;
  dissipation = 0.0;
  biomass_gibbs = 0.0;
  e_donor = -1;

  ks_coeff = nullptr;
  sub_gibbs = nullptr;
  cata_coeff = nullptr;
  anab_coeff = nullptr;
  decay_coeff = nullptr;
}


EnergyFileReader::~EnergyFileReader()
{
  delete [] line;
  delete [] keyword;
  delete [] style;
  delete [] buffer;
  memory->sfree(arg);

  if (ks_coeff != nullptr) memory->destroy(ks_coeff);
  if (sub_gibbs != nullptr) memory->destroy(sub_gibbs);
  if (cata_coeff != nullptr) memory->destroy(cata_coeff);
  if (anab_coeff != nullptr) memory->destroy(anab_coeff);
  if (decay_coeff != nullptr) memory->destroy(decay_coeff);
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::read_file(char *group_id)
{
  MPI_Barrier(world);
  int firstpass = 1;

  while (1) {
    // open file on proc 0
    if (me == 0) {
      if (firstpass)
        utils::logmesg(lmp, "Reading energy data file ...\n");
      open(filename);
    } else fp = nullptr;

    // read header info
    header(firstpass);

    // customize for new sections
    // read rest of file in free format
    while (strlen(keyword)) {

      if (strcmp(keyword, "Uptake Rate") == 0) {
        if (firstpass) uptake_rate(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Yield") == 0) {
        if (firstpass) calc_yield(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Electron Donor") == 0) {
        if (firstpass) electron_donor(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Ks Coeffs") == 0) {
        if (firstpass) ks_coeffs(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Decay Rate") == 0) {
        if (firstpass) decay_rate(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Maintenance Rate") == 0) {
        if (firstpass) maintain_rate(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Catabolic Coeffs") == 0) {
        if (firstpass) cata_coeffs(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Anabolic Coeffs") == 0) {
        if (firstpass) anab_coeffs(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Decay Coeffs") == 0) {
        if (firstpass) decay_coeffs(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Substrate Gibbs Energy") == 0) {
        if (firstpass) substrate_energy();
        else skip_lines(grid->nsubs);

      } else if (strcmp(keyword, "Biomass Gibbs Energy") == 0) {
        if (firstpass) biomass_energy(group_id);
        else skip_lines(ngroups);

      } else if (strcmp(keyword, "Dissipation Energy") == 0) {
        if (firstpass) dissipation_energy(group_id);
        else skip_lines(ngroups);

      } else
        error->all(FLERR,fmt::format("Unknown identifier in data file: {}",
                                          keyword));
      parse_keyword(0);
    }

    // close file
    if (me == 0) {
      if (compressed) pclose(fp);
      else fclose(fp);
      fp = nullptr;
    }

    // done if this was 2nd pass
    if (!firstpass) break;
    firstpass = 0;

    MPI_Barrier(world);
  }
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::uptake_rate(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    int pass = set_uptake(buf, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_uptake(char *str, char *group_id)
{
  int pass = 0;
  char *name;
  double uptake_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &uptake_one);

  if (n != 2) error->all(FLERR, "Invalid uptake rate line in data file");

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    uptake = uptake_one;
    if (uptake < 0) error->all(FLERR, "Invalid uptake rate value");
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::calc_yield(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    int pass = set_yield(buf, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_yield(char *str, char *group_id)
{
  int pass = 0;
  char *name;
  double yield_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &yield_one);

  if (n != 2) error->all(FLERR, "Invalid yield line in data file");

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    max_yield = yield_one;
    if (max_yield < 0) error->all(FLERR, "Invalid yield value");
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::electron_donor(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in Electron Donor section");
    int pass = set_edonor(arg, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_edonor(char **arg, char *group_id)
{
  int pass = 0;
  char *name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    int isub = grid->find(arg[1]);
    if (isub == -1)
      error->all(FLERR,"Invalid substrate name in Electron Donor section");
    e_donor = isub;
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::decay_rate(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    int pass = set_decay(buf, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_decay(char *str, char *group_id)
{
  int pass = 0;
  char *name;
  double decay_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &decay_one);

  if (n != 2) error->all(FLERR, "Invalid decay rate line in data file");

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    decay = decay_one;
    if (decay < 0) error->all(FLERR, "Invalid uptake rate value");
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::maintain_rate(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    int pass = set_maintain(buf, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_maintain(char *str, char *group_id)
{
  int pass = 0;
  char *name;
  double maintain_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &maintain_one);

  if (n != 2) error->all(FLERR, "Invalid maintenance rate line in data file");

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    maintain = maintain_one;
    if (decay < 0) error->all(FLERR, "Invalid maintenance rate value");
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::ks_coeffs(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  ks_coeff = memory->create(ks_coeff, grid->nsubs, "energy_file_reader:ks");

  for (int i = 0; i < grid->nsubs; i++) {
    ks_coeff[i] = 0.0;
  }

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in Ks Coeffs section");
    int pass = set_ks(narg, arg, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_ks(int narg, char **arg, char *group_id)
{
  int pass = 0;
  char *name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    int iarg = 1;
    while (iarg < narg) {
      int isub = grid->find(arg[iarg]);
      if (isub == -1)
    	error->all(FLERR,"Invalid substrate name in Ks Coeffs section");
      ks_coeff[isub] = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::cata_coeffs(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  cata_coeff = memory->create(cata_coeff, grid->nsubs, "energy_file_reader:cata_coeff");

  for (int i = 0; i < grid->nsubs; i++) {
    cata_coeff[i] = 0.0;
  }

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in Catabolic Coeffs section");
    int pass = set_cata(narg, arg, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_cata(int narg, char **arg, char *group_id)
{
  int pass = 0;
  char *name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    int iarg = 1;
    while (iarg < narg) {
      int isub = grid->find(arg[iarg]);
      if (isub == -1)
	    error->all(FLERR,"Invalid substrate name in Catabolic Coeffs section");
      cata_coeff[isub] = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::anab_coeffs(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  anab_coeff = memory->create(anab_coeff, grid->nsubs, "energy_file_reader:anab_coeff");

  for (int i = 0; i < grid->nsubs; i++) {
      anab_coeff[i] = 0.0;
  }

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in Anabolic Coeffs section");
    int pass = set_anab(narg, arg, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_anab(int narg, char **arg, char *group_id)
{
  int pass = 0;
  char *name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    int iarg = 1;
    while (iarg < narg) {
      int isub = grid->find(arg[iarg]);
      if (isub == -1)
	    error->all(FLERR,"Invalid substrate name in Anabolic Coeffs section");
      anab_coeff[isub] = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::decay_coeffs(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  decay_coeff = memory->create(decay_coeff, grid->nsubs, "energy_file_reader:decay_coeff");

  for (int i = 0; i < grid->nsubs; i++) {
    decay_coeff[i] = 0.0;
  }

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    parse_coeffs(buf);
    if (narg == 0)
      error->all(FLERR,"Unexpected empty line in Decay Coeffs section");
    int pass = set_decay_ceoffs(narg, arg, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_decay_ceoffs(int narg, char **arg, char *group_id)
{
  int pass = 0;
  char *name;
  int len = strlen(arg[0]) + 1;
  name = new char[len];
  strcpy(name,arg[0]);

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    int iarg = 1;
    while (iarg < narg) {
      int isub = grid->find(arg[iarg]);
      if (isub == -1)
	    error->all(FLERR,"Invalid substrate name in Decay Coeffs section");
      decay_coeff[isub] = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::substrate_energy()
{
  char *next;
  char *buf = new char[grid->nsubs * MAXLINE];

  sub_gibbs = memory->create(sub_gibbs, grid->nsubs, "energy_file_reader:sub_gibbs");

  for (int i = 0; i < grid->nsubs; i++) {
    sub_gibbs[i] = 0.0;
  }

  int eof = comm->read_lines_from_file(fp, grid->nsubs, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < grid->nsubs; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    set_substrate_energy(buf, i);
    buf = next + 1;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::set_substrate_energy(char *str, int i)
{
  char *name;
  double sub_gibbs_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &sub_gibbs_one);

  if (n != 2) error->all(FLERR, "Invalid substrate Gibbs energy line in data file");

  int isub = grid->find(name);
  if (isub == -1)
    error->all(FLERR,"Invalid substrate name in Substrate Gibbs Energy section");
  sub_gibbs[isub] = sub_gibbs_one;

  delete[] name;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::biomass_energy(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    int pass = set_biomass_energy(buf, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_biomass_energy(char *str, char *group_id)
{
  int pass = 0;
  char *name;
  double biomass_gibbs_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &biomass_gibbs_one);

  if (n != 2) error->all(FLERR, "Invalid biomass Gibbs energy line in data file");

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    biomass_gibbs = biomass_gibbs_one;
  }

  delete[] name;

  return pass;
}

/* ---------------------------------------------------------------------- */

void EnergyFileReader::dissipation_energy(char *group_id)
{
  char *next;
  char *buf = new char[ngroups * MAXLINE];

  int eof = comm->read_lines_from_file(fp, ngroups, MAXLINE, buf);
  if (eof) error->all(FLERR, "Unexpected end of data file");

  char *original = buf;
  for (int i = 0; i < ngroups; i++) {
    next = strchr(buf, '\n');
    *next = '\0';
    int pass = set_dissipation(buf, group_id);
    buf = next + 1;
    if (pass) break;
  }
  delete[] original;
}

/* ---------------------------------------------------------------------- */

int EnergyFileReader::set_dissipation(char *str, char *group_id)
{
  int pass = 0;
  char *name;
  double dissipation_one;
  int len = strlen(str) + 1;
  name = new char[len];

  int n = sscanf(str, "%s %lg", name, &dissipation_one);

  if (n != 2) error->all(FLERR, "Invalid dissipation energy line in data file");

  if (strcmp(group_id, name) == 0) {
    pass = 1;
    dissipation = dissipation_one;
    if (dissipation < 0) error->all(FLERR, "Invalid dissipation energy value");
  }

  delete[] name;

  return pass;
}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
------------------------------------------------------------------------- */

void EnergyFileReader::open(char *file)
{
  compressed = 0;
  char *suffix = file + strlen(file) - 3;
  if (suffix > file && strcmp(suffix, ".gz") == 0) compressed = 1;
  if (!compressed) fp = fopen(file, "r");
  else {
#ifdef LAMMPS_GZIP
    std::string gunzip = fmt::format("gzip -c -d {}",file);
#ifdef _WIN32
    fp = _popen(gunzip.c_str(),"rb");
#else
    fp = popen(gunzip.c_str(),"r");
#endif

#else
    error->one(FLERR, "Cannot open gzipped file: " + utils::getsyserror());
#endif
  }

  if (fp == nullptr)
    error->one(FLERR, fmt::format("Cannot open file {}: {}",
                                  file, utils::getsyserror()));
}

/* ----------------------------------------------------------------------
   read free-format header of data file
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
------------------------------------------------------------------------- */

void EnergyFileReader::header(int firstpass)
{
  int n;
  char *ptr;

  // customize for new sections

  const char *section_keywords[NSECTIONS] =
          {"Uptake Rate", "Yield", "Growth Rate", "Electron Donor",
           "Ks Coeffs", "Decay Rate", "Maintenance Rate",
           "Catabolic Coeffs", "Anabolic Coeffs",
           "Decay Coeffs", "Gibbs Energy", "Dissipation Energy"
          };

  // skip 1st line of file

  if (me == 0) {
    char *eof = fgets(line,MAXLINE,fp);
    if (eof == nullptr) error->one(FLERR,"Unexpected end of data file");
  }

  while (1) {

    // read a line and bcast length

    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == nullptr) n = 0;
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

    if ((ptr = strchr(line, '#'))) *ptr = '\0';
    if (strspn(line, " \t\n\r") == strlen(line)) continue;

    int rv;

    if (utils::strmatch(line, "^\\s*\\d+\\s+groups\\s")) {
      rv = sscanf(line, "%d", &ngroups);
      if (rv != 1)
        error->all(FLERR, "Could not parse 'groups' line "
                          "in data file header");
    } else break;
  }

  parse_keyword(1);
  for (n = 0; n < NSECTIONS; n++)
    if (strcmp(keyword, section_keywords[n]) == 0) break;
  if (n == NSECTIONS)
    error->all(FLERR, fmt::format("Unknown identifier in data file: {}", keyword));
}

/* ----------------------------------------------------------------------
   parse a line of coeffs into words, storing them in narg,arg
   trim anything from '#' onward
   word strings remain in line, are not copied
------------------------------------------------------------------------- */

void EnergyFileReader::parse_coeffs(char *line)
{
  char *ptr;
  if ((ptr = strchr(line,'#'))) *ptr = '\0';

  narg = 0;
  char *word = strtok(line," \t\n\r\f");

  while (word) {
    if (narg == maxarg) {
      maxarg += DELTA;
      arg = (char **)
        memory->srealloc(arg,maxarg*sizeof(char *),"read_data:arg");
    }
    arg[narg++] = word;
    word = strtok(nullptr," \t\n\r\f");
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads N lines from file
   could be skipping Natoms lines, so use bigints
------------------------------------------------------------------------- */

void EnergyFileReader::skip_lines(bigint n)
{
  if (me) return;
  if (n <= 0) return;
  char *eof = nullptr;
  for (bigint i = 0; i < n; i++) eof = fgets(line,MAXLINE,fp);
  if (eof == nullptr) error->one(FLERR,"Unexpected end of data file");
}

/* ----------------------------------------------------------------------
   grab next keyword
   read lines until one is non-blank
   keyword is all text on line w/out leading & trailing white space
   optional style can be appended after comment char '#'
   read one additional line (assumed blank)
   if any read hits EOF, set keyword to empty
   if first = 1, line variable holds non-blank line that ended header
------------------------------------------------------------------------- */

void EnergyFileReader::parse_keyword(int first)
{
  int eof = 0;
  int done = 0;

  // proc 0 reads upto non-blank line plus 1 following line
  // eof is set to 1 if any read hits end-of-file

  if (me == 0) {
    if (!first) {
      if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
    }
    while (eof == 0 && done == 0) {
      int blank = strspn(line," \t\n\r");
      if ((blank == (int)strlen(line)) || (line[blank] == '#')) {
        if (fgets(line,MAXLINE,fp) == nullptr) eof = 1;
      } else done = 1;
    }
    if (fgets(buffer,MAXLINE,fp) == nullptr) {
      eof = 1;
      buffer[0] = '\0';
    }
  }

  // if eof, set keyword empty and return

  MPI_Bcast(&eof, 1, MPI_INT, 0, world);
  if (eof) {
    keyword[0] = '\0';
    return;
  }

  // bcast keyword line to all procs

  int n;
  if (me == 0) n = strlen(line) + 1;
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  MPI_Bcast(line, n, MPI_CHAR, 0, world);

  // store optional "style" following comment char '#' after keyword

  char *ptr;
  if ((ptr = strchr(line, '#'))) {
    *ptr++ = '\0';
    while (*ptr == ' ' || *ptr == '\t') ptr++;
    int stop = strlen(ptr) - 1;
    while (ptr[stop] == ' ' || ptr[stop] == '\t'
           || ptr[stop] == '\n' || ptr[stop] == '\r')
      stop--;
    ptr[stop + 1] = '\0';
    strcpy(style, ptr);
  } else style[0] = '\0';

  // copy non-whitespace portion of line into keyword

  int start = strspn(line, " \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r')
    stop--;
  line[stop + 1] = '\0';
  strcpy(keyword, &line[start]);
}

