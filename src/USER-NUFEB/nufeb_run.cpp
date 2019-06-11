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

#include <cstring>
#include "nufeb_run.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "input.h"
#include "variable.h"
#include "compute_pressure.h"

// NUFEB specific

#include "grid.h"
#include "comm_grid.h"
#include "fix_density.h"
#include "fix_diffusion_reaction.h"
#include "fix_monod.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NufebRun::NufebRun(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg)
{
  diffdt = 1.0;
  difftol = 1.0;

  pairdt = 1.0;
  pairtol = 1.0;
  
  nfix_monod = 0;
  nfix_diff = 0;

  fix_density = NULL;
  fix_monod = NULL;
  fix_diffusion = NULL;
  comp_pressure = NULL;
  
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "diffdt") == 0) {
      diffdt = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "difftol") == 0) {
      difftol = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pairdt") == 0) {
      pairdt = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pairtol") == 0) {
      pairtol = force->numeric(FLERR, arg[iarg+1]);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal run_style nufeb command");
    }
  }
}

/* ---------------------------------------------------------------------- */

NufebRun::~NufebRun()
{
  delete [] fix_monod;
  delete [] fix_diffusion;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void NufebRun::init()
{
  // this is required because many places check for verlet style
  delete [] update->integrate_style;
  update->integrate_style = new char[7];
  strcpy(update->integrate_style, "verlet\0");

  // create fix nufeb/density
  char **fixarg = new char*[3];
  fixarg[0] = (char *)"nufeb_density";
  fixarg[1] = (char *)"all";
  fixarg[2] = (char *)"nufeb/density";
  modify->add_fix(3, fixarg, 1);
  delete [] fixarg;
  fix_density = (FixDensity *)modify->fix[modify->nfix-1];

  // allocate space for storing nufeb/monod fixes
  fix_monod = new FixMonod*[modify->nfix];
  
  // find all nufeb/monod fixes
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style, "nufeb/monod")) {
      fix_monod[nfix_monod++] = (FixMonod *)modify->fix[i];
    }
  }

  // allocate space for storing nufeb/diffusion_reaction fixes
  fix_diffusion = new FixDiffusionReaction*[modify->nfix];
  
  // find all nufeb/diffusion_reaction fixes
  for (int i = 0; i < modify->nfix; i++) {
    if (strstr(modify->fix[i]->style, "nufeb/diffusion_reaction")) {
      fix_diffusion[nfix_diff++] = (FixDiffusionReaction *)modify->fix[i];
    }
  }

  // create compute volume
  char **volarg = new char*[3];
  volarg[0] = (char *)"nufeb_volume";
  volarg[1] = (char *)"all";
  volarg[2] = (char *)"nufeb/volume";
  modify->add_compute(3, volarg, 1);
  delete [] volarg;
  comp_pressure = (ComputePressure *)modify->compute[modify->ncompute-1];

  // create compute volume variable
  char **vararg = new char*[3];
  vararg[0] = (char *)"nufeb_volume";
  vararg[1] = (char *)"equal";
  vararg[2] = (char *)"c_nufeb_volume";
  input->variable->set(3, vararg);
  delete [] vararg;
  
  // create compute pressure
  char **pressarg = new char*[8];
  pressarg[0] = (char *)"nufeb_pressure";
  pressarg[1] = (char *)"all";
  pressarg[2] = (char *)"pressure";
  pressarg[3] = (char *)"NULL";
  pressarg[4] = (char *)"pair";
  pressarg[5] = (char *)"fix";
  pressarg[6] = (char *)"vol";
  pressarg[7] = (char *)"v_nufeb_volume";
  modify->add_compute(8, pressarg, 1);
  delete [] pressarg;
  comp_pressure = (ComputePressure *)modify->compute[modify->ncompute-1];

  Integrate::init();

  // warn if no fixes

  if (modify->nfix == 0 && comm->me == 0)
    error->warning(FLERR,"No fixes defined, atoms won't move");

  // virial_style:
  // 1 if computed explicitly by pair->compute via sum over pair interactions
  // 2 if computed implicitly by pair->virial_fdotr_compute via sum over ghosts

  if (force->newton_pair) virial_style = 2;
  else virial_style = 1;

  // setup lists of computes for global and per-atom PE and pressure

  ev_setup();

  // detect if fix omp is present for clearing force arrays

  int ifix = modify->find_fix("package_omp");
  if (ifix >= 0) external_force_clear = 1;

  // set flags for arrays to clear in force_clear()

  torqueflag = extraflag = 0;
  if (atom->torque_flag) torqueflag = 1;
  if (atom->avec->forceclearflag) extraflag = 1;

  // check if simulation box is orthogonal
  if (domain->triclinic)
    error->all(FLERR, "Triclinic simulation box not supported in nufeb run style");
  triclinic = 0;
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void NufebRun::setup(int flag)
{
  if (comm->me == 0 && screen) {
    fprintf(screen,"Setting up NUFEB run ...\n");
    if (flag) {
      fprintf(screen,"  Unit style    : %s\n",update->unit_style);
      fprintf(screen,"  Current step  : " BIGINT_FORMAT "\n",update->ntimestep);
      fprintf(screen,"  Time step     : %g\n",update->dt);
      timer->print_timeout(screen);
    }
  }

  if (lmp->kokkos)
    error->all(FLERR,"KOKKOS package requires run_style nufeb/kk");

  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atom->setup();
  modify->setup_pre_exchange();
  if (triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  if (atom->sortfreq > 0) atom->sort();
  comm->borders();
  if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  domain->image_check();
  domain->box_too_small_check();
  modify->setup_pre_neighbor();
  neighbor->build(1);
  modify->setup_post_neighbor();
  neighbor->ncalls = 0;

  // NUFEB specific

  grid->setup();
  comm_grid->setup();
  
  // compute all forces

  force->setup();
  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    else force->kspace->compute_dummy(eflag,vflag);
  }

  modify->setup_pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  output->setup(flag);
  update->setupflag = 0;

  // NUFEB specific

  // compute density
  fix_density->post_integrate();
  
  // run diffusion until it reaches steady state
  int niter = diffusion();
  if (comm->me == 0)
    fprintf(screen, "Initial diffusion reaction converged in %d steps\n", niter);

  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->reaction_flag = 0;
    fix_monod[i]->growth_flag = 1;
  }

  for (int i = 0; i < nfix_diff; i++) {
    fix_diffusion[i]->compute_flag = 0;
  }
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void NufebRun::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    modify->setup_pre_exchange();
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();
    neighbor->build(1);
    modify->setup_post_neighbor();
    neighbor->ncalls = 0;
  }

  // compute all forces

  ev_set(update->ntimestep);
  force_clear();
  modify->setup_pre_force(vflag);

  if (pair_compute_flag) force->pair->compute(eflag,vflag);
  else if (force->pair) force->pair->compute_dummy(eflag,vflag);

  if (atom->molecular) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) {
    force->kspace->setup();
    if (kspace_compute_flag) force->kspace->compute(eflag,vflag);
    else force->kspace->compute_dummy(eflag,vflag);
  }

  modify->setup_pre_reverse(eflag,vflag);
  if (force->newton) comm->reverse_comm();

  modify->setup(vflag);
  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run for N steps
------------------------------------------------------------------------- */

void NufebRun::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_post_neighbor = modify->n_post_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {
    if (timer->check_timeout(i)) {
      update->nsteps = i;
      break;
    }

    ntimestep = ++update->ntimestep;

    // needs to come before ev_set
    comp_pressure->addstep(ntimestep);

    ev_set(ntimestep);
    
    growth();

    // store current dt
    double dt = update->dt;
    update->dt = pairdt;
    reset_dt();

    // disable density computation
    fix_density->compute_flag = 0;
    
    // disable compute of monod growth fixes
    for (int i = 0; i < nfix_monod; i++) {
      fix_monod[i]->compute_flag = 0;
    }

    // disable compute of diffusion fixes
    for (int i = 0; i < nfix_diff; i++) {
      fix_diffusion[i]->compute_flag = 0;
    }

    int niter = 0;
    do {
      // initial time integration

      timer->stamp();
      modify->initial_integrate(vflag);
      if (n_post_integrate) modify->post_integrate();
      timer->stamp(Timer::MODIFY);

      // regular communication vs neighbor list rebuild

      nflag = neighbor->decide();

      if (nflag == 0) {
    	timer->stamp();
    	comm->forward_comm();
    	timer->stamp(Timer::COMM);
      } else {
    	if (n_pre_exchange) {
    	  timer->stamp();
    	  modify->pre_exchange();
    	  timer->stamp(Timer::MODIFY);
    	}
    	if (triclinic) domain->x2lamda(atom->nlocal);
    	domain->pbc();
    	if (domain->box_change) {
    	  domain->reset_box();
    	  comm->setup();
    	  if (neighbor->style) neighbor->setup_bins();
    	}
    	timer->stamp();
    	comm->exchange();
    	if (sortflag && ntimestep >= atom->nextsort) atom->sort();
    	comm->borders();
    	if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    	timer->stamp(Timer::COMM);
    	if (n_pre_neighbor) {
    	  modify->pre_neighbor();
    	  timer->stamp(Timer::MODIFY);
    	}
    	neighbor->build(1);
    	timer->stamp(Timer::NEIGH);
    	if (n_post_neighbor) {
    	  modify->post_neighbor();
    	  timer->stamp(Timer::MODIFY);
    	}
      }

      // force computations
      // important for pair to come before bonded contributions
      // since some bonded potentials tally pairwise energy/virial
      // and Pair:ev_tally() needs to be called before any tallying

      force_clear();

      timer->stamp();

      if (n_pre_force) {
    	modify->pre_force(vflag);
    	timer->stamp(Timer::MODIFY);
      }

      if (pair_compute_flag) {
    	force->pair->compute(eflag,vflag);
    	timer->stamp(Timer::PAIR);
      }
    
      if (atom->molecular) {
    	if (force->bond) force->bond->compute(eflag,vflag);
    	if (force->angle) force->angle->compute(eflag,vflag);
    	if (force->dihedral) force->dihedral->compute(eflag,vflag);
    	if (force->improper) force->improper->compute(eflag,vflag);
    	timer->stamp(Timer::BOND);
      }

      if (kspace_compute_flag) {
    	force->kspace->compute(eflag,vflag);
    	timer->stamp(Timer::KSPACE);
      }

      if (n_pre_reverse) {
    	modify->pre_reverse(eflag,vflag);
    	timer->stamp(Timer::MODIFY);
      }

      // reverse communication of forces

      if (force->newton) {
    	comm->reverse_comm();
    	timer->stamp(Timer::COMM);
      }

      // force modifications, final time integration, diagnostics

      if (n_post_force) modify->post_force(vflag);
      modify->final_integrate();
      if (n_end_of_step) modify->end_of_step();
      timer->stamp(Timer::MODIFY);
    
      // double press = comp_pressure->compute_scalar();
      // fprintf(screen, "press:%e\n", press);
      ++niter;
    } while(fabs(comp_pressure->compute_scalar()) > pairtol);
    if (comm->me == 0) fprintf(screen, "pair interaction: %d steps\n", niter);

    // update densities
    fix_density->compute_flag = 1;
    fix_density->post_integrate();

    // run diffusion until it reaches steady state
    update->dt = diffdt; 
    reset_dt();
    niter = diffusion();
    // fprintf(screen, "diffusion: %d steps\n", niter);

    // restore original dt
    update->dt = dt;
    reset_dt();
    
    // all output

    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }
}

/* ---------------------------------------------------------------------- */

void NufebRun::reset_dt()
{
  if (force->pair) force->pair->reset_dt();
  for (int i = 0; i < modify->nfix; i++) modify->fix[i]->reset_dt();
}

/* ---------------------------------------------------------------------- */

void NufebRun::cleanup()
{
  modify->post_run();
  domain->box_too_small_check();
  update->update_time();
}

/* ----------------------------------------------------------------------
   clear force on own & ghost atoms
   clear other arrays as needed
------------------------------------------------------------------------- */

void NufebRun::force_clear()
{
  size_t nbytes;

  if (external_force_clear) return;

  // clear force on all particles
  // if either newton flag is set, also include ghosts
  // when using threads always clear all forces.

  int nlocal = atom->nlocal;

  if (neighbor->includegroup == 0) {
    nbytes = sizeof(double) * nlocal;
    if (force->newton) nbytes += sizeof(double) * atom->nghost;

    if (nbytes) {
      memset(&atom->f[0][0],0,3*nbytes);
      if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
      if (extraflag) atom->avec->force_clear(0,nbytes);
    }

  // neighbor includegroup flag is set
  // clear force only on initial nfirst particles
  // if either newton flag is set, also include ghosts

  } else {
    nbytes = sizeof(double) * atom->nfirst;

    if (nbytes) {
      memset(&atom->f[0][0],0,3*nbytes);
      if (torqueflag) memset(&atom->torque[0][0],0,3*nbytes);
      if (extraflag) atom->avec->force_clear(0,nbytes);
    }

    if (force->newton) {
      nbytes = sizeof(double) * atom->nghost;

      if (nbytes) {
        memset(&atom->f[nlocal][0],0,3*nbytes);
        if (torqueflag) memset(&atom->torque[nlocal][0],0,3*nbytes);
        if (extraflag) atom->avec->force_clear(nlocal,nbytes);
      }
    }
  }
}

void NufebRun::growth()
{
  // setup monod growth fixes for growth
  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->compute_flag = 1;
    fix_monod[i]->reaction_flag = 0;
    fix_monod[i]->growth_flag = 1;
  }

  // grow atoms
  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->post_integrate();
  }
}

int NufebRun::diffusion()
{
  // setup monod growth fixes for diffusion
  for (int i = 0; i < nfix_monod; i++) {
    fix_monod[i]->compute_flag = 1;
    fix_monod[i]->reaction_flag = 1;
    fix_monod[i]->growth_flag = 0;
  }

  // make sure diffusion fixes are enabled
  for (int i = 0; i < nfix_diff; i++) {
    fix_diffusion[i]->compute_flag = 1;
  }

  // store current dt before changing it
  double dt = update->dt;
  update->dt = diffdt;
  
  int niter = 0;
  bool flag;
  do {
    timer->stamp();
    comm_grid->forward_comm();
    timer->stamp(Timer::COMM);

    flag = true;
    for (int i = 0; i < nfix_diff; i++) {
      fix_diffusion[i]->initial_integrate(0);
    }
    for (int i = 0; i < nfix_monod; i++) {
      fix_monod[i]->post_integrate();
    }
    for (int i = 0; i < nfix_diff; i++) {
      fix_diffusion[i]->final_integrate();
      double res = fix_diffusion[i]->compute_scalar();
      // fprintf(screen, "%e ", res);
      flag &= res < difftol;
    }
    // fprintf(screen, "\n");
    timer->stamp(Timer::MODIFY);
    ++niter;
  } while (!flag);

  // restore previous dt
  update->dt = dt;

  return niter;
}
