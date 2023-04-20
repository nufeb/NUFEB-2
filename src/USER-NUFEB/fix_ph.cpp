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

#include <math.h>
#include <float.h>

#include "atom.h"
#include "error.h"
#include "fix_ph.h"
#include "energy_file_reader.h"
#include "grid.h"
#include "memory.h"
#include "modify.h"
#include "grid_masks.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define GAS_CONST 0.0083144   // ideal gas constant (KJ/mol.K)

/* ---------------------------------------------------------------------- */

FixPH::FixPH(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6)
    error->all(FLERR, "Illegal fix nufeb/ph command");

  temp = 298.15;
  iph = 0.0;

  sstc_gibbs = nullptr;
  ncharges = nullptr;
  form_id = nullptr;

  keq = nullptr;
  act_all = nullptr;

  buff_flag = 0;
  ina = -1;
  icl = -1;
  ih = -1;

  // read energy file
  EnergyFileReader * reader = nullptr;
  reader = new EnergyFileReader(lmp, arg[3]);
  reader->read_file(arg[1], arg[2]);

  sstc_gibbs = reader->sstc_gibbs;
  ncharges = reader->ncharges;
  form_id = reader->form_id;

  if (sstc_gibbs == nullptr || ncharges == nullptr || form_id == nullptr)
    error->all(FLERR, "Insufficient parameter sections in data file");

  iph = utils::numeric(FLERR, arg[4], true, lmp);
  if (iph < 0 || iph > 14)
    error->all(FLERR, "Illegal fix nufeb/ph command: ph");

  ih = grid->find(arg[5]);
  if (ih < 0)
    error->all(FLERR, "Can't find substrate name: H");

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temperature") == 0) {
      temp = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "buffer") == 0) {
      buff_flag = 1;
      ina = grid->find(arg[iarg + 1]);
      if (ina < 0)
        error->all(FLERR, "Can't find substrate name: na");

      icl = grid->find(arg[iarg + 2]);
      if (icl < 0)
        error->all(FLERR, "Can't find substrate name: cl");

      phlo = utils::numeric(FLERR, arg[iarg + 3], true, lmp);
      phhi = utils::numeric(FLERR, arg[iarg + 4], true, lmp);
      if (phlo < 0 || phhi > 14 || phlo > phhi)
        error->all(FLERR, "Illegal fix nufeb/ph command: phlo, phhi");
      iarg += 5;
    } else {
      error->all(FLERR, "Illegal fix nufeb/ph command");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixPH::~FixPH()
{
  memory->destroy(keq);
  memory->destroy(act_all);
  memory->destroy(sh);
}

/* ---------------------------------------------------------------------- */

int FixPH::setmask()
{
  int mask = 0;
  mask |= CHEMISTRY_NUFEB;
  mask |= REACTOR_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPH::init()
{
  int ncells = grid->ncells;
  int nsubs = grid->nsubs;

  keq = memory->create(keq, nsubs, 4, "nufeb/ph:keq");
  act_all = memory->create(act_all, nsubs, 5, ncells, "nufeb/act_all");
  sh = memory->create(sh, ncells, "nufeb/ph:sh");

  ph = memory->create(ph, ncells, "nufeb/ph:ph");
  act = memory->create(act, nsubs, ncells, "nufeb/ph:act");
  grid->ph = ph;
  grid->act = act;

  init_keq();
  compute_activity(0, ncells, iph);
}

/* ---------------------------------------------------------------------- */

void FixPH::chemistry_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixPH::reactor_nufeb()
{
  if (buff_flag) {
    to_mol();
    buffer_ph();
    to_kg();
  }
}

/* ---------------------------------------------------------------------- */

void FixPH::compute()
{
  to_mol();
  compute_ph(0, grid->ncells);
  to_kg();
}

/* ---------------------------------------------------------------------- */

void FixPH::init_keq() {
  // water Kj/mol
  double dgo_h2o = -237.18;
  double rt = temp * GAS_CONST;

  for (int i = 0; i < grid->nsubs; i++) {
    for (int j = 0; j < 4; j++) {
      int k = j+1;
      if (sstc_gibbs[i][j] != DBL_MAX && sstc_gibbs[i][k] != DBL_MAX) {
        if (j == 0) {
          // hydration equilibrium constants
          keq[i][0] = exp((dgo_h2o + sstc_gibbs[i][0] - sstc_gibbs[i][1]) / -rt);
        } else {
          // acid-base equilibrium constants
          keq[i][j] = exp((sstc_gibbs[i][k] - sstc_gibbs[i][j]) / -rt);
        }
      } else {
        keq[i][j] = 0.0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
 compute activities of five substrate protonation forms
 ------------------------------------------------------------------------- */

void FixPH::compute_activity(int first, int last, double iph) {
  int nsubs = grid->nsubs;
  int ncells = grid->ncells;
  double **conc = grid->conc;

  double *theta = memory->create(theta,nsubs,"ph:theta");
  double gsh = pow(10, -iph);
  double gsh2 = gsh * gsh;
  double gsh3 = gsh * gsh2;

  for (int i = 0; i < nsubs; i++) {
    theta[i] = (1 + keq[i][0]) * gsh3 + keq[i][1] * gsh2 + keq[i][2] * keq[i][3] * gsh
               + keq[i][3] * keq[i][2] * keq[i][1];
    if (theta[i] == 0) {
      lmp->error->all(FLERR, "theta returns a zero value");
    }
    double tmp[5];
    tmp[0] = keq[i][0] * gsh3 / theta[i];
    tmp[1] = gsh3 / theta[i];
    tmp[2] = gsh2 * keq[i][1] / theta[i];
    tmp[3] = gsh * keq[i][1] * keq[i][2] / theta[i];
    tmp[4] = keq[i][1] * keq[i][2] * keq[i][3] / theta[i];

    for (int j = first; j < last; j++) {
      if (grid->mask[j] & BLAYER_MASK) continue;
        sh[j] = gsh;
      // not hydrated form acitivity
      act_all[i][0][j] = conc[i][j] * tmp[0];
      // fully protonated form activity
      act_all[i][1][j] = conc[i][j] * tmp[1];
      // 1st deprotonated form activity
      act_all[i][2][j] = conc[i][j] * tmp[2];
      // 2nd deprotonated form activity
      act_all[i][3][j] = conc[i][j] * tmp[3];
      // 3rd deprotonated form activity
      act_all[i][4][j] = conc[i][j] * tmp[4];

      if (i == ih) {
        act_all[i][1][j] = gsh;
      }
      int form = form_id[i];
      grid->act[i][j] = act_all[i][form][j];
      grid->ph[j] = -log10(sh[j]);
    }
  }
  memory->destroy(theta);
}

/* ---------------------------------------------------------------------- */

inline double sum_activity(double ***act_all, double **keq, double **conc, int **ncharges,
                           double theta, double *gsh, int i, int j) {
  // not hydrated form acitivity
  act_all[i][0][j] = keq[i][0] * conc[i][j] * gsh[2] / theta;
  // fully protonated form activity
  act_all[i][1][j] = conc[i][j] * gsh[2] / theta;
  // 1st deprotonated form activity
  act_all[i][2][j] = conc[i][j] * gsh[1] * keq[i][1] / theta;
  // 2nd deprotonated form activity
  act_all[i][3][j] = conc[i][j] * gsh[0] * keq[i][1] * keq[i][2] / theta;
  // 3rd deprotonated form activity
  act_all[i][4][j] = conc[i][j] * keq[i][1] * keq[i][2] * keq[i][3] / theta;

  double tmp[5];
  tmp[0] = ncharges[i][0] * act_all[i][0][j];
  tmp[1] = ncharges[i][1] * act_all[i][1][j];
  tmp[2] = ncharges[i][2] * act_all[i][2][j];
  tmp[3] = ncharges[i][3] * act_all[i][3][j];
  tmp[4] = ncharges[i][4] * act_all[i][4][j];
  return tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
}


/* ---------------------------------------------------------------------- */

inline void set_gsh(double *gsh, double value) {
  gsh[0] = value;
  gsh[1] = gsh[0] * gsh[0];
  gsh[2] = gsh[1] * gsh[0];
}

/* ---------------------------------------------------------------------- */

void FixPH::compute_ph(int first, int last) {
  double tol = 5e-15;
  int max_iter = 100;

  int nsubs = grid->nsubs;
  int ncells = grid->ncells;
  double **conc = grid->conc;

  double *fa = memory->create(fa, ncells, "nufeb/ph:fa");
  double *fb = memory->create(fb, ncells, "nufeb/ph:fb");
  double *f = memory->create(f, ncells, "nufeb/ph:f");
  double *df = memory->create(df, ncells, "nufeb/ph:df");

  double a = 1e-14;
  double b = 1;

  for (int j = first; j < last; j++) {
    if (grid->mask[j] & BLAYER_MASK) continue;
    fa[j] = a;
    fb[j] = b;
  }

  double gsh[3];
  set_gsh(gsh, a);
  for (int i = 0; i < nsubs; i++) {
    double theta = (1 + keq[i][0]) * gsh[2] + keq[i][1] * gsh[1] + keq[i][2] * keq[i][1] * gsh[0]
                   + keq[i][3] * keq[i][2] * keq[i][1];
    if (theta <= 0) {
      lmp->error->all(FLERR, "theta returns a zero value");
    }
    for (int j = first; j < last; j++) {
      if (grid->mask[j] & BLAYER_MASK) continue;
      fa[j] += sum_activity(act_all, keq, conc, ncharges, theta, gsh, i, j);
    }
  }

  set_gsh (gsh, b);
  for (int i = 0; i < nsubs; i++) {
    double theta = (1 + keq[i][0]) * gsh[2] + keq[i][1] * gsh[1] + keq[i][2] * keq[i][1] * gsh[0]
                   + keq[i][3] * keq[i][2] * keq[i][1];
    if (theta <= 0) {
      lmp->error->all(FLERR, "theta returns a zero value");
    }
    for (int j = first; j < last; j++) {
      if (grid->mask[j] & BLAYER_MASK) continue;
      fb[j] += sum_activity(act_all, keq, conc, ncharges, theta, gsh, i, j);
    }
  }

  bool wrong = false;
  for (int j = first; j < last; j++) {
    if (grid->mask[j] & BLAYER_MASK) continue;
    if (fa[j] * fb[j] > 0)
      wrong = true;
  }
  if (wrong)
    lmp->error->all(FLERR, "The sum of charges returns a wrong value");

  // Newton-Raphson method
  int iter = 1;
  while (iter <= max_iter) {
    for (int j = first; j < last; j++) {
      if (grid->mask[j] & BLAYER_MASK) continue;
      f[j] = sh[j];
      df[j] = 1;
    }

    for (int i = 0; i < nsubs; i++) {
      for (int j = first; j < last; j++) {
        if (grid->mask[j] & BLAYER_MASK) continue;
        set_gsh (gsh, sh[j]);
        double theta = (1 + keq[i][0]) * gsh[2] + keq[i][1] * gsh[1] + keq[i][2] * keq[i][1] * gsh[0]
                       + keq[i][3] * keq[i][2] * keq[i][1];
        f[j] += sum_activity(act_all, keq, conc, ncharges, theta, gsh, i, j);

        double dtheta = theta * theta;
        double aux = 3 * gsh[1] * (keq[i][0] + 1) + 2 * gsh[0] * keq[i][1] + keq[i][1] * keq[i][2];
        double tmp[5];
        tmp[0] = ncharges[i][0] * ((3 * gsh[1] * keq[i][0] * conc[i][j]) / theta - (keq[i][0] * conc[i][j] * gsh[2] * aux) / dtheta);
        tmp[1] = ncharges[i][1] * ((3 * gsh[1] * conc[i][j]) / theta - (conc[i][j] * gsh[2] * aux) / dtheta);
        tmp[2] = ncharges[i][2] * ((2 * gsh[0] * keq[i][1] * conc[i][j]) / theta - (keq[i][1] * conc[i][j] * gsh[1] * aux) / dtheta);
        tmp[3] = ncharges[i][3] * ((keq[i][1] * keq[i][2] * conc[i][j]) / theta - (keq[i][1] * keq[i][2] * conc[i][j] * gsh[0] * aux) / dtheta);
        tmp[4] = ncharges[i][4] * (-(keq[i][1] * keq[i][2] * keq[i][3] * conc[i][j] * aux) / dtheta);
        df[j] += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4];
      }
    }

    // Check for convergence
    bool flag = true;
    for (int j = first; j < last; j++) {
      if (grid->mask[j] & BLAYER_MASK) continue;
      if (fabs(f[j]) >= tol)
        flag = false;
    }
    if (flag) break;

    // Compute next value
    for (int j = first; j < last; j++) {
      if (grid->mask[j] & BLAYER_MASK) continue;
      if (fabs(f[j]) >= tol) {
        double d = f[j] / df[j];
        // Prevent sh below 1e-14. That can happen because sometimes the Newton
        // method overshoots to a negative gsh value, due to a small derivative
        // value.
        if (d >= sh[j] - 1e-14)
          d = sh[j] / 2;
        sh[j] -= d;
      }
    }
    iter++;
  }

  for (int j = first; j < last; j++) {
    if (grid->mask[j] & BLAYER_MASK) continue;
    for (int i = 0; i < nsubs; i++) {
      int form = form_id[i];
      act[i][j] = act_all[i][form][j];

      if (i == ih) {
        act_all[ih][1][j] =sh[j];
        act[i][j] = sh[j];
      }
    }
    grid->ph[j] = -log10(sh[j]);
  }

  memory->destroy(fa);
  memory->destroy(fb);
  memory->destroy(f);
  memory->destroy(df);
}

/* ----------------------------------------------------------------------
 buffer ph if the value is not in defined range
 ------------------------------------------------------------------------- */

void  FixPH::buffer_ph() {
  int nsubs = grid->nsubs;
  int ncells = grid->ncells;

  double prv_sh, eva_sh, ph_unbuffer;
  double prv_act_all[5], prv_act;
  int eva = -1;
  int ext = -1;

  for (int i = ncells-1; i >= 0; i--) {
    if (grid->mask[i] & GRID_MASK) {
      if (ext < 0) {
        ext = i;
        continue;
      }
      eva = i;
      break;
    }
  }
  // domain with single grid
  if (eva < 0) {
    eva = ext;
    ext++;
  }

  prv_sh = sh[eva];
  // evaluate with dynamic ph
  compute_ph(eva, ext);
  eva_sh = sh[eva];
  ph_unbuffer = grid->ph[eva];
  sh[eva] = prv_sh;

  if (ph_unbuffer < phlo || ph_unbuffer > phhi) {
    double minus = 0;
    double plus = 0;
    compute_activity(eva, ext, iph);

    for (int i = 0; i < nsubs; i++) {
      for (int j = 0; j < 5; j++) {
        double diff, act;
        act = act_all[i][j][eva];
        diff = act * ncharges[i][j];

        if (diff > 0) plus += diff;
        else if (diff < 0) minus -= diff;
      }
      sh[eva] = prv_sh;
    }
    // convert concentrations from mol/L to kg/m3
    grid->bulk[ina] += minus;
    grid->bulk[icl] += plus + eva_sh;
  }
}

/* ---------------------------------------------------------------------- */

void FixPH::to_mol()
{
  int nsubs = grid->nsubs;
  int ncells = grid->ncells;
  double *mw = grid->mw;
  double **conc = grid->conc;

  for (int i = 0; i < nsubs; i++) {
    for (int j = 0; j < ncells; j++) {
      conc[i][j] /= mw[i];;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPH::to_kg()
{
  int nsubs = grid->nsubs;
  int ncells = grid->ncells;
  double *mw = grid->mw;
  double **conc = grid->conc;
  double **act = grid->act;

  for (int i = 0; i < nsubs; i++) {
    for (int j = 0; j < ncells; j++) {
      conc[i][j] *= mw[i];
      act[i][j] *= mw[i];
    }
  }
}