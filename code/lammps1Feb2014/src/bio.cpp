/*
 * ibm.cpp
 *
 *  Created on: 8 Nov 2016
 *      Author: bowen
 */

#include <atom.h>
#include <atom_vec_bio.h>
#include <bio.h>
#include <error.h>
#include <force.h>
#include <memory.h>
#include <string.h>
#include <cctype>
#include <cstdio>

using namespace LAMMPS_NS;

BIO::BIO(LAMMPS *lmp) : Pointers(lmp)
{
  //type
  yield = NULL;
  growth = NULL;
  ks = NULL;
  typeName = NULL;
  anabCoeff = NULL;
  catCoeff = NULL;
  typeGCoeff = NULL;
  dissipation = NULL;
  tgflag = NULL;

  //nutrient
  nnus = 0;
  iniS = NULL;
  nuName = NULL;
  diffCoeff = NULL;
  nuType = NULL;
  nuGCoeff = NULL;
  ngflag = NULL;
}

/* ---------------------------------------------------------------------- */

BIO::~BIO()
{
  memory->destroy(yield);
  memory->destroy(growth);
  memory->destroy(ks);
  memory->destroy(anabCoeff);
  memory->destroy(catCoeff);
  memory->destroy(typeGCoeff);
  memory->destroy(dissipation);
  memory->destroy(tgflag);

  memory->destroy(diffCoeff);
  memory->destroy(iniS);
  memory->destroy(nuType);
  memory->destroy(nuGCoeff);
  memory->destroy(ngflag);

  for (int i = 0; i < atom->ntypes+1; i++) {
    delete [] typeName[i];
  }

  for (int i = 0; i < nnus+1; i++) {
    delete [] nuName[i];
  }

  memory->sfree(typeName);
  memory->sfree(nuName);

}

void BIO::data_nutrients(int narg, char **arg)
{
  //printf("narg = %i, nNu = %i\n", narg, nNutrients);
  if (narg != 10) error->all(FLERR,"Incorrect args for nutrient definitions");

  int id = force->numeric(FLERR,arg[0]);
  double scell = force->numeric(FLERR,arg[3]);
  double xbcm = force->numeric(FLERR,arg[4]);
  double xbcp = force->numeric(FLERR,arg[5]);
  double ybcm = force->numeric(FLERR,arg[6]);
  double ybcp = force->numeric(FLERR,arg[7]);
  double zbcm = force->numeric(FLERR,arg[8]);
  double zbcp = force->numeric(FLERR,arg[9]);

  //nutrient name
  char *name;
  int n = strlen(arg[1]) + 1;
  name = new char[n];
  strcpy(name,arg[1]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(name[i]) && name[i] != '_')
      error->all(FLERR,"Nutrient name must be "
                 "alphanumeric or underscore characters");

  for (int i = 0; i < nnus+1; i++)
    if ((nuName[i] != NULL) && (strcmp(nuName[i], name) == 0)
        && (i != id)){
      error->one(FLERR,"Repeat nutrient names");
    }

  if (nuName[id] == NULL) {
    nuName[id] = new char[n];
  } else if (strcmp(nuName[id], name) != 0){
    error->one(FLERR,"Incompatible nutrient names");
  }

  strcpy(nuName[id],name);

  delete [] name;

  //nutrient type
  int m = strlen(arg[2]);
  if (m != 1) error->all(FLERR,"Nutrient type must be a single char, "
      "l = liq, g = gas");
  char type = arg[2][0];
  if (type == 'l') nuType[id] = 0;
  else if (type == 'g') nuType[id] = 1;
  else error->all(FLERR,"Undefined nutrient type, "
      "l = liq, g = gas");

  if (iniS == NULL) error->all(FLERR,"Cannot set nutrient concentration for this nutrient style");
  iniS[id][0] = scell;
  iniS[id][1] = xbcm;
  iniS[id][2] = xbcp;
  iniS[id][3] = ybcm;
  iniS[id][4] = ybcp;
  iniS[id][5] = zbcm;
  iniS[id][6] = zbcp;
}

/* ----------------------------------------------------------------------
   set growth values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_growth(const char *str)
{
  if (growth == NULL) error->all(FLERR,"Cannot set growth for this atom style");

  char* typeName;
  double growth_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&growth_one);

  if (n != 2) error->all(FLERR,"Invalid growth line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for growth set");

  growth[itype] = growth_one;

  if (growth[itype] < 0.0) error->all(FLERR,"Invalid growth value");

  AtomVecBio* avec = (AtomVecBio *) atom->style_match("bio");
  //set growth rate for each atom
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->type[i] == itype)
      avec->atom_mu[i] = growth[itype];
  }
}
/* ----------------------------------------------------------------------
   set ks values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_ks(const char *str)
{
  if (ks == NULL) error->all(FLERR,"Cannot set ks value for this atom style");

  char* typeName;
  double ks_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&ks_one);
  if (n != 2) error->all(FLERR,"Invalid ks line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;
  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for ks set");

  ks[itype] = ks_one;
  //mass_setflag[itype] = 1;

  if (ks[itype] < 0.0) error->all(FLERR,"Invalid ks value");
}

/* ----------------------------------------------------------------------
   set yield values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_yield(const char *str)
{
  if (yield == NULL) error->all(FLERR,"Cannot set yield for this atom style");

  char* typeName;
  double yield_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&yield_one);
  if (n != 2) error->all(FLERR,"Invalid set_yield line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_yield set");

  yield[itype] = yield_one;
  //mass_setflag[itype] = 1;

  if (yield[itype] < 0.0) error->all(FLERR,"Invalid set_yield value");
}

/* ----------------------------------------------------------------------
   set diffusion values for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_diffusion(const char *str)
{
  if (diffCoeff == NULL) error->all(FLERR,"Cannot set diffCoeff for this atom style");

  char* nuName;
  double diffu_one;
  int len = strlen(str) + 1;
  nuName = new char[len];

  int n = sscanf(str,"%s %lg",nuName,&diffu_one);
  if (n != 2) error->all(FLERR,"Invalid diffCoeff line in data file");

  int inu = find_nuID(nuName);
  delete [] nuName;

  if (inu < 1 || inu > nnus)
    error->all(FLERR,"Invalid nutrient for diffusion coefficient set");

  diffCoeff[inu] = diffu_one;
  //mass_setflag[itype] = 1;

  if (diffCoeff[inu] < 0.0) error->all(FLERR,"Invalid diffCoeff value");
}


/* ----------------------------------------------------------------------
   set dissipation values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_dissipation(const char *str)
{
  if (dissipation == NULL) error->all(FLERR,"Cannot set dissipation for this atom style");

  char* typeName;
  double diss_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&diss_one);
  if (n != 2) error->all(FLERR,"Invalid dissipation line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for dissipation set");

  dissipation[itype] = diss_one;
  //mass_setflag[itype] = 1;

  if (dissipation[itype] < 0.0) error->all(FLERR,"Invalid dissipation value");
}

/* ----------------------------------------------------------------------
   set catabolism coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_catCoeff(int narg, char **arg)
{
  if (catCoeff == NULL) error->all(FLERR,"Cannot set catCoeff for this atom style");
  if (narg != nnus+1) error->all(FLERR,"Invalid catCoeff line in data file");

  char* typeName;
  int len = strlen(arg[0]) + 1;
  typeName = new char[len];
  strcpy(typeName,arg[0]);

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for catabolism coefficient set");

  for(int i = 1; i < nnus+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    catCoeff[itype][i] = value;
  }
}

/* ----------------------------------------------------------------------
   set anabolism coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_anabCoeff(int narg, char **arg)
{
  if (anabCoeff == NULL) error->all(FLERR,"Cannot set anabCoeff for this atom style");
  if (narg != nnus+1) error->all(FLERR,"Invalid anabCoeff line in data file");

  char* typeName;
  int len = strlen(arg[0]) + 1;
  typeName = new char[len];
  strcpy(typeName,arg[0]);

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for anabolism coefficient set");

  for(int i = 1; i < nnus+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    anabCoeff[itype][i] = value;
  }
}

/* ----------------------------------------------------------------------
   set energy coefficient for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_nuGCoeff(int narg, char **arg)
{
  if (nuGCoeff == NULL) error->all(FLERR,"Cannot set energy coeff for this nutrient");
  if (narg != 7) error->all(FLERR,"Invalid nuGCOeff line in data file");

  char* nuName;
  int len = strlen(arg[0]) + 1;
  nuName = new char[len];
  strcpy(nuName,arg[0]);

  int inu = find_nuID(nuName);
  delete [] nuName;

  if (inu < 1 || inu > nnus)
    error->all(FLERR,"Invalid nutrient for nuG coefficient set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "inf") == 0) {
      nuGCoeff[inu][i] = 10001;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      nuGCoeff[inu][i] = value;
    }
  }

  int flag = atoi(arg[6]);
  if ((flag > 0) && (flag < 6) && (nuGCoeff[inu][flag-1] < 1e4))
    ngflag[inu] = flag - 1;
  else error->all(FLERR,"Invalid nutrient energy flag");
}

/* ----------------------------------------------------------------------
   set energy coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_typeGCoeff(int narg, char **arg)
{
  if (typeGCoeff == NULL) error->all(FLERR,"Cannot set energy coeff for this type");
  if (narg != 7) error->all(FLERR,"Invalid typeGCoeff line in data file");

  char* typeName;
  int len = strlen(arg[0]) + 1;
  typeName = new char[len];
  strcpy(typeName,arg[0]);

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for typeG coefficient set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "inf") == 0) {
      typeGCoeff[itype][i] = 10001;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      typeGCoeff[itype][i] = value;
    }
  }

  int flag = atoi(arg[6]);
  if ((flag > 0) && (flag < 6) && (typeGCoeff[itype][flag-1] < 1e4))
    tgflag[itype] = flag - 1;
  else error->all(FLERR,"Invalid type energy flag");
}

int BIO::find_typeID(char *name) {

  for (int i = 1; i < atom->ntypes+1; i++){
    if (typeName[i] && strcmp(typeName[i],name) == 0) {
      return i;
    }
  }
  return -1;
}

int BIO::find_nuID(char *name) {

  for (int i = 1; i < nnus+1; i++)
    if (nuName[i] && strcmp(nuName[i],name) == 0) {
      return i;
    }

  return -1;
}


