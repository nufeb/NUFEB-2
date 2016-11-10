/*
 * ibm.h
 *
 *  Created on: 8 Nov 2016
 *      Author: bowen
 */

#ifndef SRC_IBM_H_
#define SRC_IBM_H_

#include "pointers.h"

namespace LAMMPS_NS {

class BIO : protected Pointers {
 public:
  //type
  char **typeName;
  double *ks, *growth, *yield;
  double **catCoeff;
  double **anabCoeff;
  double **typeGCoeff;
  double *dissipation;

  //nutrient
  int nnus;
  char **nuName;
  int *nuType;                //nutrient types 0 = liq, 1 = gas
  double **iniS;              //inlet nutrient concentrations
  double *diffCoeff;
  double **nuGCoeff;

  BIO(class LAMMPS *);
  ~BIO();

  void data_nutrients(int, char **);
  void set_growth(const char *);
  void set_ks(const char *);
  void set_yield(const char *);
  void set_diffusion(const char *);
  void set_catCoeff(int, char **);
  void set_anabCoeff(int, char **);
  void set_nuGCoeff(int, char **);
  void set_typeGCoeff(int, char **);
  void set_dissipation(const char *);

  int find_typeID(char *name);
  int find_nuID(char *name);

};

}

#endif /* SRC_IBM_H_ */
