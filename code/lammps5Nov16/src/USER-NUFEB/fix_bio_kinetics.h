/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

#ifdef FIX_CLASS

FixStyle(kinetics,FixKinetics)

#else

#ifndef SRC_FIX_KINETICS_H
#define SRC_FIX_KINETICS_H

#include "atom.h"
#include "bio.h"
#include "fix.h"
#include "decomp_grid.h"

namespace LAMMPS_NS {

class FixKinetics : public Fix, public DecompGrid<FixKinetics> {
  friend class DecompGrid<FixKinetics>;
  friend class FixKineticsEnergy;
  friend class FixKineticsMonod;
  friend class FixKineticsThermo;
  friend class FixKineticsDiffusion;
  friend class FixKineticsPH;
  friend class FixImmigration;
  friend class DumpBio;
  friend class FixKineticsBalance;

 public:
  FixKinetics(class LAMMPS *, int, char **);
  ~FixKinetics();
  int setmask();
  virtual void pre_force(int);
  void init();
  int modify_param(int, char **);
  void migrate();

  char **var;
  int *ivar;

  int nx, ny, nz, bnz;             // number of grids in x y z axis
  int bgrids;                      // # of non-boundary grids
  int ngrids;                      // # of grids
  double iph;                      // initial ph

  double **nuS;                    // nutrient concentration [nutrient][grid]
  double **nuR;                    // nutrient consumption [nutrient][grid]
  double *nuBS;                    // concentration in boundary layer [nutrient]
  double **fV;                     // velocity field [velo][grid]
  double **qGas;                   // gas chemicals [nutrient][grid]
  double **gYield;                 // inverse yield [type][grid]
  double ***activity;              // activities of chemical species [nutrient][5 charges][grid]
  double gvol, rg;           // gas volume and gas transfer constant
  double temp, rth;                // uiversal gas constant (thermodynamics) and temperature
  double **DRGCat;                 // Gibbs free energy of catabolism [type][grid]
  double **DRGAn;                  // Gibbs free energy of anabolism [type][grid]
  double **kEq;                    // equilibrium constants [nutrient][4]
  double *sh;
  int *nuConv;
  double diffT;                    // diffusion timestep
  double blayer;
  double zhi,bzhi,zlo, xlo, xhi, ylo, yhi;
  double stepz, stepx, stepy;
  int gflag;                      // microbe growth flag 1 = update biomass; 0 = solve reaction only, growth is negligible
  int demflag;                     // flag for DEM run
  double maxheight;
  int nout;
  int niter;                        // # of iterations

  int subn[3];                     // number of grids in x y axis for this proc
  int subnlo[3],subnhi[3];         // cell index of the subdomain lower and upper bound for each axis
  double sublo[3],subhi[3];        // subdomain lower and upper bound trimmed to the grid

  Grid<double, 3> grid;
  Subgrid<double, 3> subgrid;

  class AtomVecBio *avec;
  BIO *bio;
  class FixKineticsDiffusion *diffusion;
  class FixKineticsEnergy *energy;
  class FixKineticsMonod *monod;
  class FixKineticsPH *ph;
  class FixKineticsThermo *thermo;
  class FixFluid *nufebFoam;

  void init_activity();
  void init_keq();
  void integration();
  void grow();
  double getMaxHeight();
  void update_bgrids();
  bool is_inside(int);
  int position(int);
  void reset_nuR();
  void reset_isConv();

  Subgrid<double, 3> get_subgrid() const { return subgrid; }
  int get_elem_per_cell() const;
  template <typename InputIterator, typename OutputIterator>
  OutputIterator pack_cells(InputIterator first, InputIterator last, OutputIterator result) {
    for (InputIterator it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnus; i++) {
	*result++ = nuS[i][*it];
	*result++ = nuR[i][*it];
      }
      if (energy) {
	for (int i = 1; i <= bio->nnus; i++) {
	  *result++ = qGas[i][*it];
	  for (int j = 0; j < 5; j++) {
	    *result++ = activity[i][j][*it];
	  }
	}
	for (int i = 1; i <= atom->ntypes; i++) {
	  *result++ = gYield[i][*it];
	  *result++ = DRGCat[i][*it];
	  *result++ = DRGAn[i][*it];
	}
      }
      if (nufebFoam) {
	for (int i = 0; i < 3; i++) {
	  *result++ = fV[i][*it];
	}
      }
    }
    return result;
  }
  template <typename InputIterator0, typename InputIterator1>
  InputIterator1 unpack_cells(InputIterator0 first, InputIterator0 last, InputIterator1 input) {
    for (InputIterator0 it = first; it != last; ++it) {
      for (int i = 1; i <= bio->nnus; i++) {
	nuS[i][*it] = *input++;
	nuR[i][*it] = *input++;
      }
      if (energy) {
	for (int i = 1; i <= bio->nnus; i++) {
          qGas[i][*it] = *input++;
	  for (int j = 0; j < 5; j++) {
	    activity[i][j][*it] = *input++;
	  }
	}
	for (int i = 1; i <= atom->ntypes; i++) {
	  gYield[i][*it] = *input++;
	  DRGCat[i][*it] = *input++;
	  DRGAn[i][*it] = *input++;
	}
      }
      if (nufebFoam) {
	for (int i = 0; i < 3; i++) {
	  fV[i][*it] = *input++;
	}
      }
    }
    return input;
  }
  void resize(const Subgrid<double, 3> &);
};

}

#endif
#endif

