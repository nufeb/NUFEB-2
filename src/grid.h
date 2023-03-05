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

#ifndef LMP_GRID_H
#define LMP_GRID_H

#include "pointers.h"
#include <map>

enum {DIRICHLET, NEUMANN, PERIODIC};

namespace LAMMPS_NS {

  class Grid : protected Pointers {
  public:
    char *grid_style;
    class GridVec *gvec;

    typedef GridVec *(*GridVecCreator)(LAMMPS *);
    typedef std::map<std::string, GridVecCreator> GridVecCreatorMap;
    GridVecCreatorMap *gvec_map;

    bool grid_exist;
    int nmax;
    int nsubs;                  // # of substrates
    char **sub_names;           // substrate names
    double cell_size;
    int box[3];                 // # of global cells in each dimension
    int extbox[3];              // # of extended cells in each dimension
    int sublo[3], subhi[3];     // sub-box bounds in grid coordinates
    int subbox[3];              // # of cells on this proc in each dimension
    int ncells;                 // total # of cells
    int periodic[3];            // flag if x, y and z boundaries are periodic

    Grid(class LAMMPS *);
    virtual ~Grid();
    void modify_params(int, char **);
    void create_gvec(const char *, int, char **, int);
    virtual class GridVec *new_gvec(const char *, int, int &);
    void init();
    void setup();
    int find(const char *);
    int cell(double *);

    int *mask;

    // nufeb/simple
    int simple_flag;
    double **conc;    // concentration

    // nufeb/chemostat
    int chemostat_flag;
    double **reac;    // reaction rate
    double **dens;    // density
    double **diff_coeff; //diffusion Coefficient
    double ***growth; // growth rate
    double *bulk;     // bulk concentration
    double *mw;       // molecular weight g/mol
    int **boundary;   // boundary conditions (-x, +x, -y, +y, -z, +z)

  private:
    template<typename T>
    static GridVec *gvec_creator(LAMMPS *);
  };

}

#endif

/* ERROR/WARNING messages:

*/
