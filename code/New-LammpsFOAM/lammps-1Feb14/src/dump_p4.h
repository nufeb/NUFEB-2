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

DumpStyle(p4,DumpP4)

#else

#ifndef LMP_DUMP_P4_H
#define LMP_DUMP_P4_H

#include "dump.h"


namespace LAMMPS_NS {

class DumpP4 : public Dump {
 public:
  DumpP4(class LAMMPS *, int, char**);
  ~DumpP4();

 private:
  void init_style();
  void write_header(bigint);
  int count();
  void pack(int *);
  void openfile();
  void write();
  void write_data(int, double *);
  void write_conts(int, double *);
  void write_walls(int, double *);

  class Fix **fix;
  int nfix;

  char *filename_p4p;    // p3c file
  FILE *fp4p;            // file to write p3c dump to
  char *filename_p4c;    // p3c file
  FILE *fp4c;            // file to write p3c dump to
  char *format4c;
  char *filename_p4w;    // p3w file
  FILE *fp4w;            // file to write p3w dump to
  char *format4w;
 
  class Compute *compute_pairs;
  int nconts;
  int size_conts;
  int nwalls;
  int size_walls;

  double *buf4c;
  int maxbuf4c;
  bigint ctotal;
  
  double *buf4w;
  int maxbuf4w;
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
