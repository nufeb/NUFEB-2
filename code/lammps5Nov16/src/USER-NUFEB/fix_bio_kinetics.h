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

#include "fix.h"

namespace LAMMPS_NS {

class FixKinetics : public Fix {
  friend class FixKineticsEnergy;
  friend class FixKineticsMonod;
  friend class FixKineticsThermo;
  friend class FixKineticsDiffusion;
  friend class FixKineticsPH;
  friend class FixImmigration;
  friend class DumpBio;

 public:
  FixKinetics(class LAMMPS *, int, char **);
  ~FixKinetics();
  int setmask();
  virtual void pre_force(int);
  void init();

  char **var;
  int *ivar;

  int nx, ny, nz, bnz;             // number of grids in x y z axis
  int bgrids;                      // # of non-boundary grids
  int ngrids;                      // # of grids
  double iph;                      // initial ph
  int nnus;                        // # of nutrients

  double **nuS;                    // nutrient concentration [nutrient][grid]
  double **nuR;                    // nutrient consumption [nutrient][grid]
  double **fV;                     // velocity field [velo][grid]
  double **qGas;                   // gas chemicals [nutrient][grid]
  double **gYield;                 // inverse yield [type][grid]
  double ***activity;              // activities of chemical species [nutrient][5 charges][grid]
  double gVol, gasTrans;           // gas volume and gas transfer constant
  double temp, rth;                // uiversal gas constant (thermodynamics) and temperature
  double **DRGCat;                 // Gibbs free energy of catabolism [type][grid]
  double **DRGAn;                  // Gibbs free energy of anabolism [type][grid]
  double **kEq;                    // equilibrium constants [nutrient][4]
  double *Sh;
  bool *nuConv;
  double diffT;                    // diffsuion timestep
  double bl;
  double zhi,bzhi,zlo, xlo, xhi, ylo, yhi;
  double stepz, stepx, stepy;
  int gflag;                      // microbe growth flag 1 = update biomass; 0 = solve reaction only, growth is negligible

  int subn[2];                     // number of grids in x y axis for this proc
  double sublo[2], subhi[2];       // subdomain size trimmed to the grid

  class AtomVecBio *avec;
  class BIO *bio;
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
  int position(int);
  void reset_nuR();
  void reset_isConv();
};

}

#endif
#endif

