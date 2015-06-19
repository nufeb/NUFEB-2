/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_nugrowth.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixNuGrowth::FixNuGrowth(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 17) error->all(FLERR,"Illegal fix growth command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix growth command");

  var = new char*[13];
  ivar = new int[13];

  int i;
  for (i = 0; i < 13; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }
}

/* ---------------------------------------------------------------------- */

FixNuGrowth::~FixNuGrowth()
{
  int i;
  for (i = 0; i < 13; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixNuGrowth::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNuGrowth::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix growth requires atom attribute diameter");

  int i;
  for (i = 0; i < 13; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix nugrowth does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix nugrowth is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void FixNuGrowth::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  change_dia();
}

/* ---------------------------------------------------------------------- */

void FixNuGrowth::change_dia()
{

  modify->clearstep_compute();

  double KsHET = input->variable->compute_equal(ivar[0]);
  double Ko2HET = input->variable->compute_equal(ivar[1]);
  double Kno2HET = input->variable->compute_equal(ivar[2]);
  double Kno3HET = input->variable->compute_equal(ivar[3]);
  double Knh4AOB = input->variable->compute_equal(ivar[4]);
  double Ko2AOB = input->variable->compute_equal(ivar[5]);
  double Kno2NOB = input->variable->compute_equal(ivar[6]);
  double Ko2NOB = input->variable->compute_equal(ivar[7]);
  double MumHET = input->variable->compute_equal(ivar[8]);
  double MumAOB = input->variable->compute_equal(ivar[9]);
  double MumNOB = input->variable->compute_equal(ivar[10]);
  double etaHET = input->variable->compute_equal(ivar[11]);
  // double bHET = input->variable->compute_equal(ivar[12]);
  // double bAOB = input->variable->compute_equal(ivar[13]);
  // double bNOB = input->variable->compute_equal(ivar[14]);
  double bEPS = input->variable->compute_equal(ivar[12]);

  double density;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *sub = atom->sub;
  double *o2 = atom->o2;
  double *nh4 = atom->nh4;
  double *no2 = atom->no2;
  double *no3 = atom->no3;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      double gHET = 0;
      double gAOB = 0;
      double gNOB = 0;
      double gEPS = 0;
      if (type[i] == 1) {
        gHET = 1;
      }
      if (type[i] == 2) {
        gAOB = 1;
      }
      if (type[i] == 3) {
        gNOB = 1;
      }
      if (type[i] == 4) {
        gEPS = 1;
      }

      double R1 = MumHET*(sub[i]/(KsHET+sub[i]))*(o2[i]/(Ko2HET+o2[i]));
      double R2 = MumAOB*(nh4[i]/(Knh4AOB+nh4[i]))*(o2[i]/(Ko2AOB+o2[i]));
      double R3 = MumNOB*(no2[i]/(Kno2NOB+no2[i]))*(o2[i]/(Ko2NOB+o2[i]));
      double R4 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no3[i]/(Kno3HET+no3[i]))*(Ko2HET/(Ko2HET+o2[i]));
      double R5 = etaHET*MumHET*(sub[i]/(KsHET+sub[i]))*(no2[i]/(Kno2HET+no2[i]))*(Ko2HET/(Ko2HET+o2[i]));
      // double R6 = bHET;
      // double R7 = bAOB;
      // double R8 = bNOB;
      double R9 = bEPS;

      double value = update->dt * (gHET*(R1+R4+R5) + gAOB*R2 + gNOB*R3 - gEPS*R9);
      
      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      double oldMass = rmass[i];
      rmass[i] = rmass[i]*(1 + (value*nevery));
      if (rmass[i] <= 0) {
        rmass[i] = oldMass;
      }
      radius[i] = pow(((6*rmass[i])/(density*MY_PI)),(1.0/3.0))*0.5;
    }

  }

  modify->addstep_compute(update->ntimestep + nevery);
}
