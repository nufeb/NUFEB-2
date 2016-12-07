/*
 * bio.h
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
  //type (species)
  char **typeName;            // type name

  double *ks;                 // half-saturation constant
  double *mu;                 // maximum species growth rate
  double *yield;              // growth yield coefficient
  double *dissipation;        // universal gas constant (thermodynamics)

  double **catCoeff;          // catabolism coefficient [type][nutrient]
  double **anabCoeff;         // anabolism coefficient [type][nutrient]
  double **typeGCoeff;        // Gibbs free energy coefficient [type][5charges]
  int *tgflag;                // Gibbs free energy flag

  //nutrient
  int nnus;                   // # of nutrients
  int *nuType;                // nutrient types 0 = liq, 1 = gas
  char **nuName;              // nutrient name

  double *diffCoeff;          // diffusion coefficient
  double **iniS;              // inlet nutrient concentrations [nutrient][1grid + 5bc]
  double **nuGCoeff;          // Gibbs free energy coefficient [nutrient][5charges]
  int *ngflag;                // Gibbs free energy flag

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
  void set_group_mask();

  int find_typeID(char *name);
  int find_nuID(char *name);

};

}

#endif /* SRC_IBM_H_ */
