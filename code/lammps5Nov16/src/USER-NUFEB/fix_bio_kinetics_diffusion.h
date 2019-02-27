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
  friend DecompGrid<FixKineticsDiffusion> ;

public:
  FixKineticsDiffusion(class LAMMPS *, int, char **);
  ~FixKineticsDiffusion();

  char **var;
  int *ivar;

  double *rmass;
  int ntypes;
  int *mask;
  int *type;

  double stepx, stepy, stepz;             // grid size

  int xbcflag, ybcflag, zbcflag;          // boundary condition flag, 0=PERIODIC-PERIODIC, 1=DIRiCH-DIRICH, 2=NEU-DIRICH, 3=NEU-NEU, 4=DIRICH-NEU
  int bulkflag;                           // 1=solve mass balance for bulk liquid
  int shearflag, dragflag, dcflag;        // flags for shear, drag(nufebfoam), and diffusion coefficent

  double srate;                           // shear rate
  double tol;                             // tolerance for convergence criteria

  double **nugrid;                        // nutrient concentration in ghost grid [nutrient][grid], unit in mol or kg/m3
  double **xgrid;                         // grid coordinate [gird][3]
  bool *ghost;                            // ghost grid flag [gird] 1=ghost gird, 0=non-ghost grid
  double **nuprev;                        // nutrient concentration in previous diffusion step
  double **grid_diff_coeff;               // diffusion coeffs at each grid

  double vol;                             // grid volume

  double q, rvol, af;                     // parameters used for dynamic bulk
  int unit;                               // concentration unit 0=mol/l; 1=kg/m3

  int nx, ny, nz;                         // # of non-ghost grids in x, y and z
  int nxx, nyy, nzz;                      // # of all grids in x, y and z
  int nxx_yy_zz;                          // total # of grids
  double diff_dt;                         // diffusion timestep
  double xlo, xhi, ylo, yhi, zlo, zhi, bzhi;
  double xbcm, xbcp, ybcm, ybcp, zbcm, zbcp; // inlet BC concentrations for each surface

  MPI_Request *requests;

  BIO *bio;
  FixKinetics *kinetics;
  AtomVecBio *avec;

  int setmask();
  void init();
  int *diffusion(int*, int, double);
  void update_nus();
  void update_grids();
  void update_diff_coeff();
  void init_grid();
  void compute_bc(double &, double *, int, double);
  void compute_bulk();
  void compute_blayer();
  void compute_flux(double, double &, double *, double, int, int);

  bool is_equal(double, double, double);
  int get_index(int);
  void migrate(const Grid<double, 3> &, const Box<int, 3> &, const Box<int, 3> &);

  int get_elem_per_cell() const;
  template<typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnu; i++) {
        *result++ = nugrid[i][*it];
      }
    }
    return result;
  }
  template<typename InputIterator0, typename InputIterator1>
  InputIterator1 unpack_cells(InputIterator0 first, InputIterator0 last, InputIterator1 input) {
    for (InputIterator0 it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnu; i++) {
        nugrid[i][*it] = *input++;
      }
    }
    return input;
  }
  void resize(const Subgrid<double, 3> &);
};

}

#endif
#endif

