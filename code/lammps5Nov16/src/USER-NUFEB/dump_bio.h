/* -*- c++ -*- ----------------------------------------------------------
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

DumpStyle(bio,DumpBio)

#else

#ifndef LMP_DUMP_BIO_H
#define LMP_DUMP_BIO_H

#include <mpi.h>
#include "dump.h"

namespace LAMMPS_NS {

class DumpBio : public Dump {
 public:
  DumpBio(LAMMPS *, int, char**);
  virtual ~DumpBio();

 protected:
  int nevery;                // dump frequency for output
  //gzFile gzFp;  // file pointer for the compressed output stream

  FILE *fp;                  // file to write dump to
  char *filename;            // user-specified file
  int nkeywords;

  int nfix;                  // # of Fix objects used by dump
  char **id_fix;             // their IDs
  class Fix **fix;           // list of ptrs to the Fix objects

  int anFlag, concFlag, avgsFlag, catFlag, phFlag, massFlag, gasFlag, yieldFlag;
  int diaFlag, dimFlag, divFlag, heightFlag, roughFlag, segFlag, ntypeFlag, bulkFlag;
  int massHeader, divHeader, typeHeader, bulkHeader, avgsHeader;

  double xlo,xhi,ylo,yhi,zlo,zhi;
  double stepx, stepy, stepz;
  int nx, ny, nz;               // # of grids in x, y and z

  char **keywords;
  class FixKinetics *kinetics;
  class BIO *bio;

  int ncompute;                  // # of Compute objects used by dump
  class ComputeNufebDiameter *cdia;
  class ComputeNufebDimension *cdim;
  class ComputeNufebDiversity *cdiv;
  class ComputeNufebHeight *cheight;
  class ComputeNufebRough *crough;
  class ComputeNufebSegregate *cseg;
  class ComputeNufebNtypes *ctype;
  class ComputeNufebBiomass *cmass;
  class ComputeNufebAvgcon *cavgs;


  void openfile();
  void write();
  void init_style();
  void write_header(bigint);
  void pack(tagint *);
  void write_data(int, double *);

  void write_concentration_data(int);
  void write_avgcon_data();
  void write_gas_data(int);
  void write_DGRCat_data(int);
  void write_DGRAn_data(int);
  void write_pH_data();
  void write_yield_data(int);
  void write_biomass_data();
  void write_ntype_data();
  void write_bulk_data();

  void write_diameter_data();
  void write_dimension_data();
  void write_diversity_data();
  void write_height_data();
  void write_rough_data();
  void write_segregate_data();

  bigint memory_usage();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump movie filename

The file produced by dump movie cannot be binary or compressed
and must be a single file for a single processor.

E: Support for writing movies not included

LAMMPS was not built with the -DLAMMPS_FFMPEG switch in the Makefile

E: Failed to open FFmpeg pipeline to file %s

The specified file cannot be opened.  Check that the path and name are
correct and writable and that the FFmpeg executable can be found and run.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
