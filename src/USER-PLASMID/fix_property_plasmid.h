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

#ifdef FIX_CLASS

FixStyle(nufeb/property/plasmid,FixPropertyPlasmid)

#else

#ifndef LMP_FIX_PROPERTY_PLASMID_H
#define LMP_FIX_PROPERTY_PLASMID_H

#include "fix_property.h"
#include "atom_vec_bacillus.h"

#include <iostream>
#include <fstream>

namespace LAMMPS_NS {

class FixPropertyPlasmid : public FixProperty {
 public:

  FixPropertyPlasmid(class LAMMPS *, int, char **);
  ~FixPropertyPlasmid();

  int setmask();
  void grow_arrays(int);
  void biology_nufeb();
  void init() {};
  void set_arrays(int) {};
  void update_arrays(int, int);
  void compute();

  void get_plasmid_coords(int, int, double *, double *);
  void get_plasmid_coords(int, int, double *);
  void relocate_plm_x(int, int);
  void delete_filament(int *, int);
  void set_plm_x(int, int, double *,double *);

  double memory_usage();
  void copy_arrays(int, int, int);
  void copy_plasmid(int, int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  int plm_max;               // maximum # of plasmid
  int plm_init;              // # of initial plasmid
  double plm_dia;            // plasmid diameter
  double **plm_x;            // plasmid position
  double **pre_x;            // previous cell location

  double **nproteins;        // number of initiator proteins

  int fila_max;              // maximum # of filaments
  int *nfilas;               // # of filaments
  int ***fila;               // filament defined as straight line
  double **tfila;            // filament duration

  int rep_flag, par_flag;

  std::ofstream myfile;

 private:
  void get_cell_boundary(double *, int, int);
  void get_quat(double *, double *, double *);
  void distance_bt_pt_line(double *, double *, double *, double &);
  //void dump();

  class AtomVecBacillus *avec;

  double dot(double x1, double x2, double x3, double y1, double y2, double y3) { return x1*y1 + x2*y2 + x3*y3; }
  double norm(double x1, double x2, double x3) { return sqrt(dot(x1, x2, x3, x1, x2, x3)); }
  double dist(double u1, double u2, double u3, double v1, double v2, double v3) { return norm(u1-v1, u2-v2, u3-v3); }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix property/atom mol when atom_style already has molecule attribute

Self-explanatory.

E: Fix property/atom cannot specify mol twice

Self-explanatory.

E: Fix property/atom q when atom_style already has charge attribute

Self-explanatory.

E: Fix property/atom cannot specify q twice

Self-explanatory.

E: Fix property/atom rmass when atom_style already has rmass attribute

UNDOCUMENTED

E: Fix property/atom cannot specify rmass twice

UNDOCUMENTED

E: Fix property/atom vector name already exists

The name for an integer or floating-point vector must be unique.

W: Fix property/atom mol or charge or rmass w/out ghost communication

UNDOCUMENTED

E: Atom style was redefined after using fix property/atom

This is not allowed.

E: Incorrect %s format in data file

A section of the data file being read by fix property/atom does
not have the correct number of values per line.

E: Too few lines in %s section of data file

Self-explanatory.

E: Invalid atom ID in %s section of data file

An atom in a section of the data file being read by fix property/atom
has an invalid atom ID that is <= 0 or > the maximum existing atom ID.

U: Fix property/atom mol or charge w/out ghost communication

A model typically needs these properties defined for ghost atoms.

*/
