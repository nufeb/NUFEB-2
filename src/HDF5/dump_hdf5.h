/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
DumpStyle(nufeb/hdf5, DumpHDF5)
#else
#ifndef LMP_DUMP_HDF5_H
#define LMP_DUMP_HDF5_H

#include "dump.h"

extern "C" {
#include <hdf5.h>
}

#include <string>
#include <vector>

namespace LAMMPS_NS {

class DumpHDF5 : public Dump {
 public:
  DumpHDF5(LAMMPS *, int, char **);
  ~DumpHDF5();

protected:
  void write();
  void init_style() {};
  void setup();
  void write_header(bigint) {}
  void pack(tagint *) {}
  void write_data(int, double *) {}
  int parse_fields(int narg, char **arg);
  void create_one_file();
  double memory_usage() {return 0;}

  hid_t create_filespace_atom(bool oneperproc);
  hid_t create_filespace_grid(bool oneperproc);

  template <class T> 
  void trim_grids(T *bufin, T *bufout);
  template <class T>
  herr_t write_atoms_scalar(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset);
  template <class T>
  herr_t write_atoms_comp(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset, int comp);
  template <class T>
  herr_t write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc);
  template <class T>
  herr_t write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int n, int index);

  std::vector<std::string> fields;

  int ncells;
  int subdims[3], substart[3], dims[3];
};
}

#endif // LMP_DUMP_HDF5_H
#endif // DUMP_CLASS
