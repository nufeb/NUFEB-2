/*
 * fix_kinetics/diffusionS.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics/diffusion,FixKineticsDiffusion)

#else

#ifndef SRC_FIX_KINETICS_DIFFUSIONS_H
#define SRC_FIX_KINETICS_DIFFUSIONS_H

#include "bio.h"
#include "fix.h"
#include "decomp_grid.h"

namespace LAMMPS_NS {

class FixKineticsDiffusion: public Fix, protected DecompGrid<FixKineticsDiffusion> {
  friend DecompGrid<FixKineticsDiffusion>;

 public:
  FixKineticsDiffusion(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion();
  int setmask();
  void init();
  int *diffusion(int*, int, double);
  void update_nuS();

  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU

 private:
  char **var;
  int *ivar;
  int nnus;                     // # of nutrients
  double stepx, stepy, stepz;

  double *rmass;
  int ntypes;                       // # of species
  int *mask;
  int *type;
  int shearflag, dragflag;

  double shearRate;
  double **iniS;                // inlet concentrations of nutrients
  double *diffCoeff;            // diffusion coefficients of nutrients
  double *mw;                   // molecular weights of nutrients
  double tol;                   // tolerance for convergence criteria for nutrient balance equation
 // int rstep;                    // steps skipped between Si+n-Si
 // int rflag;
  double **nuR;
  double **nuS;
  double *diffD;

  double bl;
  double q, rvol, af;
  bool reactor;
  int unit;                     // nutrient unit 0 = mol/l; 1 = kg/m3
  double *nuBS;

  int nx, ny, nz;               // # of non-ghost grids in x, y and z
  int nX, nY, nZ;               // # of all grids in x, y and z
  int nXYZ;                   // total # of grids
  double diffT;
  double xlo, xhi, ylo, yhi, zlo, zhi, bzhi;
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; // inlet BC concentrations for each surface

  double **xGrid;                     // grid coordinate
  double **nuGrid;                    // nutrient concentration in ghost mesh [grid][nutrient]
  double **nuPrev;
  bool *ghost;

  MPI_Request *requests;

  Grid<double, 3> grid;
  Box<int, 3> subgrid;

  BIO *bio;
  class FixKinetics *kinetics;
  class AtomVecBio *avec;

  void update_grids();
  void compute_bc(double &, double *, int, double);
  void compute_bulk(int);
  void compute_bl();
  void compute_flux(double, double &, double *, double, int);
  bool isEuqal(double, double, double);
  int get_index(int);
  Grid<double, 3> get_grid() const;
  Box<int, 3> get_subgrid() const;
  std::array<bool, 3> get_periodic_boundary() const;
  template <typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnus; i++) {
	*result++ = nuGrid[i][*it];
      }
    }
    return result;
  }
  template <typename InputIterator0, typename InputIterator1>
  InputIterator1 unpack_cells(InputIterator0 first, InputIterator0 last, InputIterator1 input) {
    for (InputIterator0 it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnus; i++) {
	nuGrid[i][*it] = *input++;
      }
    }
    return input;
  }
  int get_cell_data_size(int n);
};

}

#endif
#endif

