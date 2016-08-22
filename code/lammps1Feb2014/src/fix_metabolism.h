/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(metabolism,FixMetabolism)

#else

#ifndef SRC_FIX_METABOLISM_H
#define SRC_FIX_METABOLISM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMetabolism : public Fix {
 public:
  FixMetabolism(class LAMMPS *, int, char **);
  ~FixMetabolism();
  int setmask();
  void init();
  void pre_force(int);

 private:
  int me;
  char *line,*keyword,*buffer;
  int compressed;
  FILE *fp;

  char **var;
  int *ivar;
  int nsubs;
  char **subName;
  int diffevery;
  int outputevery;
  double *xCell;
  double *yCell;
  double *zCell;
  double *cellVol;
  bool *ghost;

  double *cellConc;
  double *bcConc;

  int **catCoeff;
  int **anabCoeff;

  double *diffCoeff;
  double diffT;

  int numCells;
  int nx, ny, nz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  int bflag; // 1 = dirichlet, 2 = neumann, 3 = mixed
  double xstep, ystep, zstep;

  void change_dia();


  void readfile(char*);
  void open(char *);
  void header();
  void allocate_arrays();
  void parse_keyword(int first);
  void skip_lines(bigint n);

  void substrate();
  void diffusion();
  void cat();
  void anab();

  void set_diffusion(const char *str);
  void set_substrate(const char *str);
  void set_cat(const char *str);
  void set_anab(const char *str);

  int find_subID(char *name);
};

}

#endif
#endif

