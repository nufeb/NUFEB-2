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
DumpStyle(grid/vtk, DumpGridVTK)
#else

#ifndef LMP_DUMP_GRID_VTK_H
#define LMP_DUMP_GRID_VTK_H

#include "dump.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include <functional>
#include <vector>

namespace LAMMPS_NS {

class DumpGridVTK : public Dump {
 public:
  DumpGridVTK(LAMMPS *, int, char **);
  ~DumpGridVTK();

 protected:
  void write();
  void init_style();
  void write_header(bigint) {}
  void pack(tagint *) {}
  void write_data(int, double *) {}
  int parse_fields(int narg, char **arg);
  double memory_usage() {return 0;}

  void pack_concentration(vtkSmartPointer<vtkImageData>);
  void pack_reaction(vtkSmartPointer<vtkImageData>);
  void pack_density(vtkSmartPointer<vtkImageData>);
  void pack_growth(vtkSmartPointer<vtkImageData>);
  void pack_tuple1(vtkSmartPointer<vtkImageData>, const char *, double *);
  void pack_tuple1(vtkSmartPointer<vtkImageData>, const char *, double **, char **, int);
  template <int>
  void pack_tuple(vtkSmartPointer<vtkImageData>, const char *, double ***, char **, int);

  std::vector<std::string> fields;
  std::vector<std::function<void(vtkSmartPointer<vtkImageData>)> > packs;
};
}

#endif
#endif
