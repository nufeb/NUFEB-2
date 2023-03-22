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
#include "atom.h"
#include "error.h"

#include "atom_vec_bacillus.h"
#include "fix_growth_energy.h"
#include "energy_file_reader.h"
#include "grid.h"
#include "memory.h"
#include "modify.h"
#include "grid_masks.h"
#include "math_const.h"
#include "update.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define MAINT_ENERGY 0.00125 // maintenance energy KJ/mol·s, derived from 4.5KJ/mol·h

/* ---------------------------------------------------------------------- */

FixGrowthEnergy::FixGrowthEnergy(LAMMPS *lmp, int narg, char **arg) :
  FixGrowth(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR, "Illegal fix nufeb/growth/energy command");

  temp = 298.15;
  gas_const = 0.0083144;
  alfa = 1.2;
  beta = 0.8;
  mw_biomass = 24.6;
  dgo_cata = 0.0;
  dgo_anab = 0.0;

  uptake = 0.0;
  decay = 0.0;
  dissipation = 0.0;
  biomass_gibbs = 0.0;
  e_donor = -1;

  gibbs_cata = nullptr;
  gibbs_anab = nullptr;
  yield = nullptr;

  sub_gibbs = nullptr;
  ks_coeff = nullptr;
  cata_coeff = nullptr;
  anab_coeff = nullptr;
  decay_coeff = nullptr;

  // read energy file
  EnergyFileReader * reader = nullptr;
  reader = new EnergyFileReader(lmp, arg[3]);
  reader->read_file(arg[1]);

  uptake = reader->uptake;
  decay = reader->decay;
  dissipation = reader->dissipation;
  biomass_gibbs = reader->biomass_gibbs;
  sub_gibbs = reader->sub_gibbs;
  ks_coeff = reader->ks_coeff;
  cata_coeff = reader->cata_coeff;
  anab_coeff = reader->anab_coeff;
  decay_coeff = reader->decay_coeff;
  e_donor = reader->e_donor;

  if (sub_gibbs == nullptr || ks_coeff == nullptr || cata_coeff == nullptr ||
      anab_coeff == nullptr || decay_coeff == nullptr)
    error->all(FLERR, "Insufficient parameter sections defined in energy file");

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "temperature") == 0) {
      temp = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if  (strcmp(arg[iarg], "gas_const") == 0) {
      gas_const = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if  (strcmp(arg[iarg], "alfa") == 0) {
      gas_const = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if  (strcmp(arg[iarg], "beta") == 0) {
      gas_const = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else if  (strcmp(arg[iarg], "mw_biomass") == 0) {
      mw_biomass = utils::numeric(FLERR, arg[iarg + 1], true, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/growth/energy command");
    }
  }

  if (alfa < beta)
    error->all(FLERR, "Illegal alfa and beta values in fix nufeb/growth/energy command");

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
}

/* ---------------------------------------------------------------------- */

FixGrowthEnergy::~FixGrowthEnergy() {
  memory->destroy(gibbs_cata);
  memory->destroy(gibbs_anab);
  memory->destroy(yield);
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::init()
{
  if (grid->mw == nullptr)
    error->all(FLERR, "fix nufeb/growth/energy requires molecular weight (mw) defined");

  dt = update->dt;

  gibbs_cata = memory->create(gibbs_cata, grid->ncells, "growth/energy:gibbs_cata");
  gibbs_anab = memory->create(gibbs_anab, grid->ncells, "growth/energy:gibbs_anab");
  yield = memory->create(yield, grid->ncells, "growth/energy:yield");

  compute_dgo();
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::biology_nufeb()
{
  if (update->ntimestep % nevery) return;
    compute_dgr();
    update_atoms();
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::chemistry_nufeb()
{
  compute_dgr();
  update_cells();
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::compute_dgo()
{
  for (int i = 0; i < grid->nsubs; i++) {
    dgo_cata += cata_coeff[i] * sub_gibbs[i];
    dgo_anab += anab_coeff[i] * sub_gibbs[i];
  }
  dgo_cata += biomass_gibbs;
  dgo_anab += biomass_gibbs;
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::compute_dgr()
{
  double **conc = grid->conc;
  double *mw = grid->mw;
  double rt = temp * gas_const;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GHOST_MASK) continue;
    // standard Gibbs free energy
    gibbs_cata[i] = dgo_cata;
    gibbs_anab[i] = dgo_anab;

    double q_cata = 1.0;
    double q_anab = 1.0;

    // compute catabolic and anabolic reaction quotient Q
    for (int j = 0; j < grid->nsubs; j++) {
      if (conc[j][i] <= 0) continue;

      // convert concentrations from kg/m3 to mol/L
      if (cata_coeff[j] > 0) {
        q_cata *= pow(conc[j][i]/mw[j], cata_coeff[j]);
      } else if (cata_coeff[j] < 0){
        q_cata /= pow(conc[j][i]/mw[j], -cata_coeff[j]);
      }

      if (anab_coeff[j] > 0) {
        q_anab *= pow(conc[j][i]/mw[j], anab_coeff[j]);
      } else if (anab_coeff[j] < 0){
        q_anab /= pow(conc[j][i]/mw[j], -anab_coeff[j]);
      }
    }

    // Gibbs free energy
    gibbs_cata[i] += rt * log(q_cata);
    gibbs_anab[i] += rt * log(q_anab);
  }
}

/* ---------------------------------------------------------------------- */

double FixGrowthEnergy::compute_monod(int cell)
{
  double monod = 1.0;
  double **conc = grid->conc;

  for (int i = 0; i < grid->nsubs; i++) {
    if (ks_coeff[i] > 0)
      monod *= conc[i][cell] / (ks_coeff[i] + conc[i][cell]);
  }

  return monod;
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::update_cells()
{
  double **reac = grid->reac;
  double **dens = grid->dens;
  double *mw = grid->mw;

  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GHOST_MASK) continue;
    double inv_yield;
    double spec_growth;
    double m_req, q_cat;
    double meta_coeff;

    yield[i] = -gibbs_cata[i] / (gibbs_anab[i] + dissipation);
    // inverse yield
    if (yield[i] != 0.0)
      inv_yield = 1 / yield[i];
    else
      inv_yield = 0.0;

    // Specific substrate uptake rate for catabolism
    // unit mole-eD/molCx·s
    q_cat = uptake * compute_monod(i);
    // specific substrate consumption required for maintenance
    // unit = mol-eD / mol-X·s
    m_req = -MAINT_ENERGY / gibbs_cata[i];
    spec_growth = (q_cat - m_req) * yield[i];

    // update substrate reaction term
    for (int j = 0; j < grid->nsubs; j++) {
      if (q_cat > alfa * m_req) {
        meta_coeff = inv_yield * cata_coeff[j] + anab_coeff[j];
        // convert unit mol-eD/mol-X to kg-eD/kg-X
        meta_coeff *= mw[e_donor] / mw_biomass;
        reac[j][i] += spec_growth * meta_coeff * dens[igroup][i];

      } else if (q_cat <= alfa * m_req && q_cat >= beta * m_req) {
        reac[j][i] += spec_growth * cata_coeff[j] * dens[igroup][i];

      } else if (q_cat < beta * m_req) {
        reac[j][i] -= spec_growth * decay_coeff[j] * dens[igroup][i];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGrowthEnergy::update_atoms()
{
  for (int i = 0; i < grid->ncells; i++) {
    if (grid->mask[i] & GHOST_MASK) continue;
    double spec_growth;
    double m_req, q_cat;
    // yield unit = mol-X/mol-eD
    yield[i] = -gibbs_cata[i] / (gibbs_anab[i] + dissipation);

    // specific substrate uptake rate for catabolism
    // unit mole-eD/mol-x·s
    q_cat = uptake * compute_monod(i);
    // specific substrate consumption required for maintenance
    // unit = mol-eD / mol-X·s
    m_req = -MAINT_ENERGY / gibbs_cata[i];
    spec_growth = (q_cat - m_req) * yield[i];

    if (q_cat > alfa * m_req) {
      grid->growth[igroup][i][0] = spec_growth;

    } else if (q_cat <= alfa * m_req && q_cat > beta * m_req) {
      grid->growth[igroup][i][0] = 0.0;

    } else if (q_cat <= beta * m_req) {
      grid->growth[igroup][i][0] = -decay * (m_req - q_cat) / m_req;
    }
  }

  if (atom->coccus_flag) {
    update_atoms_coccus();
  } else {
    update_atoms_bacillus(avec);
  }
}

