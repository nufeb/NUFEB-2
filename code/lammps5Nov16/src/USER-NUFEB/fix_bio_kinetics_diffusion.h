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
class AtomVecBio;
class BIO;
class FixKinetics;

class FixKineticsDiffusion: public Fix, public DecompGrid<FixKineticsDiffusion> {
  friend DecompGrid<FixKineticsDiffusion>;

 public:
  FixKineticsDiffusion(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion();

  char **var;
  int *ivar;

  double stepx, stepy, stepz;

  int xbcflag, ybcflag, zbcflag;             // 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  int bulkflag;                           // 1=solve mass balance for bulk liquid

  double *rmass;
  int ntypes;                       // # of species
  int *mask;
  int *type;
  int shearflag, dragflag;

  double srate;                 // shear rate
  double tol;                   // tolerance for convergence criteria for nutrient balance equation

  double **nuGrid;                    // nutrient concentration in ghost mesh [nutrient][grid]
  double **xGrid;                     // grid coordinate [gird][3]
  bool *ghost;

  double *diffD;

  double bl;
  double q, rvol, af;
  bool reactor;
  int unit;                     // nutrient unit 0 = mol/l; 1 = kg/m3

  int nx, ny, nz;               // # of non-ghost grids in x, y and z
  int nX, nY, nZ;               // # of all grids in x, y and z
  int nXYZ;                   // total # of grids
  double diffT;
  double xlo, xhi, ylo, yhi, zlo, zhi, bzhi;
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; // inlet BC concentrations for each surface

  double **nuPrev;

  MPI_Request *requests;

  BIO *bio;
  FixKinetics *kinetics;
  AtomVecBio *avec;

  int setmask();
  void init();
  int *diffusion(int*, int, double);
  void update_nuS();
  void update_grids();
  void init_grid();
  void compute_bc(double &, double *, int, double);
  void compute_bulk();
  void compute_bl();
  void compute_flux(double, double &, double *, double, int, int);
  bool isEuqal(double, double, double);
  int get_index(int);
  void migrate(const Grid<double, 3> &, const Box<int, 3> &, const Box<int, 3> &);

  int get_elem_per_cell() const;
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
  void resize(const Subgrid<double, 3> &);
};

}

#endif
#endif

