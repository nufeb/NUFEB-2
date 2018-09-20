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
  //type (microbial species)
  char **tname;            // type name

  double **ks;                // half-saturation constant [type][nutrient]
  double *q;                  // specific consumption rate
  double *mu;                 // specific growth rate
  double *yield;              // growth yield coefficient
  double *dissipation;        // universal gas constant (thermodynamics)
  double *maintain;           // maintenance [type]
  double *decay;              // decay rate [type]
  double *edoner;

  double **cata_coeff;          // catabolism coefficient [type][nutrient]
  double **anab_coeff;         // anabolism coefficient [type][nutrient]
  double **decay_coeff;        // decay coefficient [type][nutrient]
  double **tgibbs_coeff;        // Gibbs free energy coefficient [type][5charges]
  int *tgibbs_flag;                // Gibbs free energy flag for type [type]
  int **tcharge;              // charge [type][5charges]

  //nutrient
  int nnu;                   // # of nutrients
  int *nustate;                // nutrient types 0 = liq, 1 = gas
  char **nuname;              // nutrient name

  double *diff_coeff;          // diffusion coefficient [nutrient]
  double *mw;                 // molecular Weights [nutrient]
  double **ini_nus;              // inlet nutrient concentrations [nutrient][1grid + 5bc]
  double **nugibbs_coeff;          // Gibbs free energy coefficient [nutrient][5charges]
  int *ngflag;                // Gibbs free energy flag for nutrients [nutrient]
  int **nucharge;                // charge [nutrient][5charges]
  double *kla;                // mass Transfer Coefficient [nutrient]

  BIO(class LAMMPS *);
  ~BIO();

  void type_grow();
  void create_type(char *);
  void data_nutrients(int, char **);
  void set_typeName(int narg, char **arg);
  void set_q(const char *);
  void set_mu(const char *);
  void set_mw(const char *);
  void set_ks(int, char **);
  void set_yield(const char *);
  void set_edoner(const char *);
  void set_maintain(const char *);
  void set_decay(const char *);
  void set_diffusion(const char *);
  void set_catCoeff(int, char **);
  void set_decayCoeff(int, char **);
  void set_anabCoeff(int, char **);
  void set_nuGCoeff(int, char **);
  void set_typeGCoeff(int, char **);
  void set_dissipation(const char *);
  void set_nuChr(int, char **);
  void set_typeChr(int, char **);
  void set_kLa(const char *);
  void set_group_mask();

  int find_typeID(char *name);
  int find_nuID(char *name);

};

}

#endif /* SRC_IBM_H_ */
