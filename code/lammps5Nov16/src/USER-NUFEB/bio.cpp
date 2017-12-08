/*
 * ibm.cpp
 *
 *  Created on: 8 Nov 2016
 *      Author: bowen
 */

#include "bio.h"

#include <string.h>
#include <cctype>
#include <cstdio>
#include <cstdlib>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "atom_vec_bio.h"

using namespace LAMMPS_NS;

BIO::BIO(LAMMPS *lmp) : Pointers(lmp)
{
  //type
  yield = NULL;
  maintain = NULL;
  decay = NULL;
  q = NULL;
  mu = NULL;
  ks = NULL;
  eD = NULL;
  typeName = NULL;
  anabCoeff = NULL;
  catCoeff = NULL;
  decayCoeff = NULL;
  typeGCoeff = NULL;
  dissipation = NULL;
  tgflag = NULL;
  typeChr = NULL;

  //nutrient
  nnus = 0;
  iniS = NULL;
  nuName = NULL;
  diffCoeff = NULL;
  nuType = NULL;
  nuGCoeff = NULL;
  ngflag = NULL;
  nuChr = NULL;
  kLa = NULL;
  mw = NULL;
}

/* ---------------------------------------------------------------------- */

BIO::~BIO()
{
  memory->destroy(yield);
  memory->destroy(maintain);
  memory->destroy(decay);
  memory->destroy(eD);
  memory->destroy(q);
  memory->destroy(mu);
  memory->destroy(ks);
  memory->destroy(anabCoeff);
  memory->destroy(catCoeff);
  memory->destroy(decayCoeff);
  memory->destroy(typeGCoeff);
  memory->destroy(dissipation);
  memory->destroy(tgflag);
  memory->destroy(typeChr);

  memory->destroy(diffCoeff);
  memory->destroy(mw);
  memory->destroy(iniS);
  memory->destroy(nuType);
  memory->destroy(nuGCoeff);
  memory->destroy(ngflag);
  memory->destroy(nuChr);
  memory->destroy(kLa);

  for (int i = 0; i < atom->ntypes+1; i++) {
    delete [] typeName[i];
  }

  for (int i = 0; i < nnus+1; i++) {
    delete [] nuName[i];
  }

  memory->sfree(typeName);
  memory->sfree(nuName);
}

void BIO::type_grow()
{
  if (yield != NULL) memory->grow(yield,atom->ntypes+1,"bio:yield");
  if (maintain != NULL) memory->grow(maintain,atom->ntypes+1,"bio:maintain");
  if (decay != NULL) memory->grow(decay,atom->ntypes+1,"bio:decay");
  if (dissipation != NULL) memory->grow(dissipation,atom->ntypes+1,"bio:dissipation");
  if (q != NULL) memory->grow(q,atom->ntypes+1,"bio:q");
  if (mu != NULL) memory->grow(mu,atom->ntypes+1,"bio:mu");
  if (ks != NULL) memory->grow(ks,atom->ntypes+1,nnus+1,"bio:ks");
  if (anabCoeff != NULL) memory->grow(anabCoeff,atom->ntypes+1,nnus+1,"bio:anabCoeff");
  if (catCoeff != NULL) memory->grow(catCoeff,atom->ntypes+1,nnus+1,"bio:catCoeff");
  if (decayCoeff != NULL) memory->grow(decayCoeff,atom->ntypes+1,nnus+1,"bio:decayCoeff");
  if (typeGCoeff != NULL) memory->grow(typeGCoeff,atom->ntypes+1,5,"bio:typeGCoeff");
  if (tgflag != NULL) memory->grow(tgflag,atom->ntypes+1,"bio:tgflag");
  if (typeChr != NULL) memory->grow(typeChr,atom->ntypes+1,5,"bio:typeChr");
  if (eD != NULL) memory->grow(eD,atom->ntypes+1,"bio:eD");
  if (typeName != NULL) typeName = (char **) memory->srealloc(typeName,(atom->ntypes+1)*sizeof(char *),"bio:typeName");
}

void BIO::create_type(char *name) {
  atom->ntypes = atom->ntypes + 1;
  type_grow();
  typeName[atom->ntypes] = name;

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

void BIO::set_q(const char *str)
{
  if (q == NULL) error->all(FLERR,"Cannot set growth rate for this atom style");

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

  q[itype] = growth_one;

  if (q[itype] < 0.0) error->all(FLERR,"Invalid growth value");

 // AtomVecBio* avec = (AtomVecBio *) atom->style_match("bio");
  //set consumption rate for each atom
//  for (int i = 0; i < atom->nlocal; i++) {
//    if (atom->type[i] == itype)
//      avec->atom_q[i] = q[itype];
//  }
}

/* ----------------------------------------------------------------------
   set growth values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_mu(const char *str)
{
  if (mu == NULL) error->all(FLERR,"Cannot set consumption rate for this atom style");

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

  mu[itype] = growth_one;

  if (mu[itype] < 0.0) error->all(FLERR,"Invalid growth value");

  AtomVecBio* avec = (AtomVecBio *) atom->style_match("bio");
}

/* ----------------------------------------------------------------------
   set ks values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_ks(int narg, char **arg)
{
  if (ks == NULL) error->all(FLERR,"Cannot set Ks for this atom style");
  if (narg != nnus+1) error->all(FLERR,"Invalid Ks line in data file");

  char* typeName;
  int len = strlen(arg[0]) + 1;
  typeName = new char[len];
  strcpy(typeName,arg[0]);

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for Ks set");

  for(int i = 1; i < nnus+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    ks[itype][i] = value;
  }
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

  if (yield_one < 0)
    lmp->error->all(FLERR,"yield cannot be zero or less than zero");

  yield[itype] = yield_one;
  //mass_setflag[itype] = 1;

  if (yield[itype] < 0.0) error->all(FLERR,"Invalid set_yield value");
}


/* ----------------------------------------------------------------------
   set yield values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_eD(const char *str)
{
  if (eD == NULL) error->all(FLERR,"Cannot set eD for this atom style");

  char* typeName;
  double eD_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&eD_one);
  if (n != 2) error->all(FLERR,"Invalid set_eD line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_eD set");

  eD[itype] = eD_one;
  //mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set maintenance values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_maintain(const char *str)
{
  if (maintain == NULL) error->all(FLERR,"Cannot set maintain for this atom style");

  char* typeName;
  double maintain_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&maintain_one);
  if (n != 2) error->all(FLERR,"Invalid set_maintain line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_maintain set");

  maintain[itype] = maintain_one;
  //mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set decay values for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_decay(const char *str)
{
  if (decay == NULL) error->all(FLERR,"Cannot set decay for this atom style");

  char* typeName;
  double decay_one;
  int len = strlen(str) + 1;
  typeName = new char[len];

  int n = sscanf(str,"%s %lg",typeName,&decay_one);
  if (n != 2) error->all(FLERR,"Invalid set_decay line in data file");

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for set_decay set");

  decay[itype] = decay_one;
  //mass_setflag[itype] = 1;
}

/* ----------------------------------------------------------------------
   set diffusion coeff for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_diffusion(const char *str)
{
  if (diffCoeff == NULL) error->all(FLERR,"Cannot set diffCoeff for this nutrient style");

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
   set molecular weights for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_mw(const char *str)
{
  if (mw == NULL) error->all(FLERR,"Cannot set molecular weights for this nutrient style");

  char* nuName;
  double mw_one;
  int len = strlen(str) + 1;
  nuName = new char[len];

  int n = sscanf(str,"%s %lg",nuName,&mw_one);
  if (n != 2) error->all(FLERR,"Invalid molecular weights line in data file");

  int inu = find_nuID(nuName);
  delete [] nuName;

  if (inu < 1 || inu > nnus)
    error->all(FLERR,"Invalid nutrient for molecular weights set");

  mw[inu] = mw_one;
  //mass_setflag[itype] = 1;

  if (mw[inu] < 0.0) error->all(FLERR,"Invalid molecular weights value");
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
   set decay coefficient for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_decayCoeff(int narg, char **arg)
{
  if (decayCoeff == NULL) error->all(FLERR,"Cannot set decayCoeff for this atom style");
  if (narg != nnus+1) error->all(FLERR,"Invalid decayCoeff line in data file");

  char* typeName;
  int len = strlen(arg[0]) + 1;
  typeName = new char[len];
  strcpy(typeName,arg[0]);

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for decay coefficient set");

  for(int i = 1; i < nnus+1; i++) {
    double value = force->numeric(FLERR,arg[i]);
    decayCoeff[itype][i] = value;
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

/* ----------------------------------------------------------------------
   set charges for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_nuChr(int narg, char **arg)
{
  if (nuChr == NULL) error->all(FLERR,"Cannot set charge for this nutrient");
  if (narg != 6) error->all(FLERR,"Invalid nutrient charge line in data file");

  char* nuName;
  int len = strlen(arg[0]) + 1;
  nuName = new char[len];
  strcpy(nuName,arg[0]);

  int inu = find_nuID(nuName);
  delete [] nuName;

  if (inu < 1 || inu > nnus)
    error->all(FLERR,"Invalid nutrient for nutrient charge set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "na") == 0) {
      nuChr[inu][i] = 0;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      nuChr[inu][i] = value;
    }
  }
}

/* ----------------------------------------------------------------------
   set charges for all types
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_typeChr(int narg, char **arg)
{
  if (typeChr == NULL) error->all(FLERR,"Cannot set charge for this type");
  if (narg != 6) error->all(FLERR,"Invalid type charge line in data file");

  char* typeName;
  int len = strlen(arg[0]) + 1;
  typeName = new char[len];
  strcpy(typeName,arg[0]);

  int itype = find_typeID(typeName);
  delete [] typeName;

  if (itype < 1 || itype > atom->ntypes)
    error->all(FLERR,"Invalid type for typeG coefficient set");

  for(int i = 0; i < 5; i++) {
    if (strcmp(arg[i+1], "na") == 0) {
      typeChr[itype][i] = 0;
    } else {
      double value = force->numeric(FLERR,arg[i+1]);
      typeChr[itype][i] = value;
    }
  }
}

/* ----------------------------------------------------------------------
   set mass transfer coeff for all nutrients
   called from reading of data file
------------------------------------------------------------------------- */

void BIO::set_kLa(const char *str)
{
  if (kLa == NULL) error->all(FLERR,"Cannot set KLa for this atom style");

  char* nuName;
  double kLa_one;
  int len = strlen(str) + 1;
  nuName = new char[len];

  int n = sscanf(str,"%s %lg",nuName,&kLa_one);
  if (n != 2) error->all(FLERR,"Invalid KLa line in data file");

  int inu = find_nuID(nuName);
  delete [] nuName;

  if (inu < 1 || inu > nnus)
    error->all(FLERR,"Invalid nutrient for KLa set");

  kLa[inu] = kLa_one;
  //mass_setflag[itype] = 1;

  if (kLa[inu] < 0.0) error->all(FLERR,"Invalid KLa value");
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
