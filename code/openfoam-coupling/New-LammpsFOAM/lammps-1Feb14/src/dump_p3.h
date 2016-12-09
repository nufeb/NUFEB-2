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

#ifdef DUMP_CLASS

DumpStyle(p3,DumpP3)

#else

#ifndef LMP_DUMP_P3_H
#define LMP_DUMP_P3_H

#include "dump.h"
#include <stdlib.h>     /* atoi */


namespace LAMMPS_NS {

class DumpP3 : public Dump {
 public:
  DumpP3(class LAMMPS *, int, char**);
  ~DumpP3();

 private:
  void init_style();
  void write_header(bigint);
  void parse_extra_variables(int,char**);
  void split_var_name(char*,char*,char*,int&);
  void add_compute(char*,int&,int&);
  void add_fix(char*,int&,int&);
  void init_extra_variables_pointer();
  int count();
  void pack(int *);
  void pack_style_p3(int *);
  void pack_style_p4(int *);
  void openfile();
  void write();
  void write_data(int, double *);
  void write_conts(int, double *);
  void write_walls(int, double *);

  int nevery;

  class Fix **fix;
  char **fix_name;
  int nfix_names;
  bool fix_all;
  int nfix;

  int n_extra_vars;
  double **var_ptr;      // pointer to extra vars
  int     *var_stride;   // stride for extra vars access
  int     *var_type;     // ATOM, COMPUTE or FIX
  int     *var_index;    // local index for arrays
  char    *header_extra; // List of variables in .p3p file
  int     n_atom_vars;
  int     n_compute_vars;
  int     n_fix_vars;
  int     *var_atom;     // enum with possible vars (omegax,fy,...)
  int     *var_compute;  // ...
  int     *var_fix;      // ...

  class Compute **vcomputes;
  char **id_vcomputes;
  int   n_vcomputes;

  class Fix **vfixs;
  char **id_vfixs;
  int   n_vfixs;

  int      p3_style;     // P3 or P4 style


  char *filename_p3p;    // p3c file
  FILE *fp3p;            // file to write p3c dump to
  char *filename_p3c;    // p3c file
  FILE *fp3c;            // file to write p3c dump to
  char *format3c;
  char *filename_p3w;    // p3w file
  FILE *fp3w;            // file to write p3w dump to
  char *format3w;
 
  bool use_p3c;
  bool use_p3w;

  class Compute *compute_pairs;
  char *compute_name;
  int nconts;
  int size_conts;
  int nwalls;
  int size_walls;

  double *buf3c;
  int maxbuf3c;
  bigint ctotal;
  
  double *buf3w;
  int maxbuf3w;
  bigint wtotal;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid dump p3 filename

Filenames used with the dump p3 style cannot be binary or cause files
to be written by each processor.

*/
