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
#include "grid.h"

namespace LAMMPS_NS {

class FixKinetics : public Fix {
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
  void borders();
  int modify_param(int, char **);

  char **var;
  int *ivar;

  int nx, ny, nz, bnz;             // number of grids in x y z axis
  int bgrids;                      // # of non-boundary grids
  int ngrids;                      // # of grids
  double iph;                      // initial ph
  int nnus;                        // # of nutrients

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

  int subn[3];                     // number of grids in x y axis for this proc
  int subnlo[3],subnhi[3];         // cell index of the subdomain lower and upper bound for each axis
  double sublo[3],subhi[3];        // subdomain lower and upper bound trimmed to the grid

  int recv_buff_size;
  int send_buff_size;
  int *recvcells;
  int *sendcells;
  int *recvbegin, *recvend;
  int *sendbegin, *sendend;

  int *cellbegin;                  // first atom index for each cell
  int *next;                       // points to the next atom index lying in the same cell

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
  void update_bgrids();
  bool is_inside(int);
  int position(int);
  void add_cells(const Grid &, const Grid &, int *, int);
  bool is_intersection_valid(const Grid &);
  void send_recv_cells(const Grid &, const Grid &, const Grid &, int &, int &);
  void reset_nuR();
  void reset_isConv();
#ifdef OUTPUT_GRID
  void output_nutrient_info(int);
  void output_bacteria_info(int);
#endif
};

}

#endif
#endif

