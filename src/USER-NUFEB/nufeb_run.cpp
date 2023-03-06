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

#ifdef _WIN32
#include <windows.h>
#include <stdint.h> // <cstdint> requires C++-11
#else
#include <sys/time.h>
#include <sys/resource.h>
#endif

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
#include "compute_ke.h"

// NUFEB specific

#include "grid.h"
#include "comm_grid.h"
#include "fix_density.h"
#include "fix_diffusion_reaction.h"
#include "compute_volume.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NufebRun::NufebRun(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg)
{
  init_diff_flag = true;
  info = true;
  ndiff = 0;
  npair = 0;

  biodt = 1.0;
  diffdt = 1.0;
  difftol = 1.0;
  diffmax = -1;
  pairdt = 1.0;
  pairtol = 1.0;
  pairmax = -1;

  nfix_diffusion = 0;

  fix_density = nullptr;
  fix_diffusion = nullptr;
  comp_pressure = nullptr;
  comp_ke = nullptr;
  comp_volume = nullptr;

  profile = nullptr;
  
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "diffdt") == 0) {
      diffdt = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "difftol") == 0) {
      difftol = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "diffmax") == 0) {
      diffmax = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pairdt") == 0) {
      pairdt = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pairtol") == 0) {
      pairtol = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pairmax") == 0) {
      pairmax = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "profile") == 0) {
      char filename[80];
      sprintf(filename, "%s_%d.log", arg[iarg+1], comm->me);
      profile = fopen(filename,"w");
      iarg += 2;
    } else if (strcmp(arg[iarg], "screen") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) info = true;
      else if (strcmp(arg[iarg+1], "no") == 0) info = false;
      else {
	error->all(FLERR, "Illegal run_style nufeb command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "initdiff") == 0) {
      if (strcmp(arg[iarg+1], "yes") == 0) init_diff_flag = true;
      else if (strcmp(arg[iarg+1], "no") == 0) init_diff_flag = false;
      else {
	error->all(FLERR, "Illegal run_style nufeb command");
      }
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal run_style nufeb command");
    }
  }
}

/* ---------------------------------------------------------------------- */

NufebRun::~NufebRun()
{
  if (profile)
    fclose(profile);
  
  delete [] fix_diffusion;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void NufebRun::init()
{
  // this is required because many places check for verlet style
  delete [] update->integrate_style;
  update->integrate_style = new char[13];
  strcpy(update->integrate_style, "verlet/nufeb\0");

  // allocate space for storing fix diffusion
  fix_diffusion = new FixDiffusionReaction*[modify->nfix];
  
  for (int i = 0; i < modify->nfix; i++) {
    // find nufeb fixes
    if (strstr(modify->fix[i]->style, "nufeb/diffusion_reaction")) {
      fix_diffusion[nfix_diffusion++] = (FixDiffusionReaction *)modify->fix[i];
    }
  }

  // create fix nufeb/density
  char **fixarg = new char*[3];
  fixarg[0] = (char *)"nufeb_density";
  fixarg[1] = (char *)"all";
  fixarg[2] = (char *)"nufeb/density";
  modify->add_fix(3, fixarg, 1);
  delete [] fixarg;
  fix_density = (FixDensity *)modify->fix[modify->nfix-1];

  // create compute volume
  char **volarg = new char*[3];
  volarg[0] = (char *)"nufeb_volume";
  volarg[1] = (char *)"all";
  volarg[2] = (char *)"nufeb/volume";
  modify->add_compute(3, volarg, 1);
  delete [] volarg;
  comp_volume = (ComputeVolume *)modify->compute[modify->ncompute-1];
  
  // create compute pressure
  char **pressarg = new char*[6];
  pressarg[0] = (char *)"nufeb_pressure";
  pressarg[1] = (char *)"all";
  pressarg[2] = (char *)"pressure";
  pressarg[3] = (char *)"NULL";
  pressarg[4] = (char *)"pair";
  pressarg[5] = (char *)"fix";
  modify->add_compute(6, pressarg, 1);
  delete [] pressarg;
  comp_pressure = (ComputePressure *)modify->compute[modify->ncompute-1];

  // create compute ke
  char **kearg = new char*[3];
  kearg[0] = (char *)"nufeb_ke";
  kearg[1] = (char *)"all";
  kearg[2] = (char *)"ke";
  modify->add_compute(3, kearg, 1);
  delete [] kearg;
  comp_ke = (ComputeKE *)modify->compute[modify->ncompute-1];

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

  // set timestep to biological dt

  biodt = update->dt;

  // compute density
  fix_density->compute();

  // run diffusion until it reaches steady state
  if (init_diff_flag) {
    int niter = module_chemsitry();
    if (comm->me == 0)
      fprintf(screen, "Initial diffusion reaction converged in %d steps\n", niter);
  } else {
    if (comm->me == 0)
      fprintf(screen, "Initial diffusion reaction convergence disabled\n");
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

  for (int i = 0; i < n; i++) {
    double step_start = get_time();

    if (timer->check_timeout(i)) {
      update->nsteps = i;
      break;
    }

    ntimestep = ++update->ntimestep;

    // needs to come before ev_set
    comp_pressure->addstep(ntimestep);

    ev_set(ntimestep);
    timer->stamp();

    // run biology module
    double t = get_time();
    module_biology();
    if (profile)
      fprintf(profile, "%d %e ", update->ntimestep, get_time()-t);
    
    // run physics module
    t = get_time();
    double press;
    press = module_physics();
    if (profile)
      fprintf(profile, "%d %e ", npair, get_time()-t);
    if (info && comm->me == 0) fprintf(screen, "pair interaction: %d steps (pressure %e N/m2)\n", npair, press);

    // reset to biological timestep
    update->dt = biodt;
    reset_dt();

    // call all fixes implementing post_physics_nufeb()
    if (modify->n_post_physics_nufeb) {
      timer->stamp();
      modify->post_physics_nufeb();
      timer->stamp(Timer::MODIFY);
    }

    // run chemistry module
    t = get_time();
    ndiff = module_chemsitry();
    if (profile)
      fprintf(profile, "%d %e ", ndiff, get_time()-t);
    if (info && comm->me == 0) fprintf(screen, "diffusion: %d steps\n", ndiff);

    // run reactor module
    t = get_time();
    module_reactor();
    if (profile)
      fprintf(profile, "%e\n", get_time()-t);

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

/* ---------------------------------------------------------------------- */

void NufebRun::module_biology()
{
  // reset to biological timestep

  update->dt = biodt;
  reset_dt();

  // call all fixes implementing biology_nufeb()

  if (modify->n_biology_nufeb) {
    timer->stamp();
    modify->biology_nufeb();
    timer->stamp(Timer::MODIFY);
  }
}

/* ---------------------------------------------------------------------- */

int NufebRun::module_chemsitry()
{
  // reset to diffusion timestep

  update->dt = diffdt;
  reset_dt();
  
//  for (int i = 0; i < nfix_diffusion; i++) {
//    fix_diffusion[i]->closed_system_initial();
//  }

  int niter = 0;
  bool conv_flag;
  bool converge[nfix_diffusion];
  for (int i = 0; i < nfix_diffusion; i++) {
    converge[i] = false;
  }

  do {
    timer->stamp();
    comm_grid->forward_comm();
    timer->stamp(Timer::COMM);

    for (int i = 0; i < nfix_diffusion; i++) {
      fix_diffusion[i]->compute_initial();
    }

    // call all fixes implementing chemistry_nufeb()
    if (modify->n_chemistry_nufeb) {
      timer->stamp();
      modify->chemistry_nufeb();
      timer->stamp(Timer::MODIFY);
    }

    conv_flag = true;
    for (int i = 0; i < nfix_diffusion; i++) {
      if (converge[i]) continue;
      fix_diffusion[i]->compute_final();
      double res = fix_diffusion[i]->compute_scalar();
      if (res < difftol) converge[i] = true;
      if (!converge[i]) conv_flag = false;
    }

    timer->stamp(Timer::MODIFY);
    ++niter;

    if (diffmax > 0 && niter >= diffmax)
      conv_flag = true;

  } while (!conv_flag);

  for (int i = 0; i < nfix_diffusion; i++) {
    fix_diffusion[i]->closed_system_scaleup(biodt);
  }

  return niter;
}

/* ---------------------------------------------------------------------- */

double NufebRun::module_physics()
{
  bigint ntimestep;
  int nflag,sortflag;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_post_neighbor = modify->n_post_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  // reset to DEM timestep

  update->dt = pairdt;
  reset_dt();

  double vol = comp_volume->compute_scalar();
  timer->stamp(Timer::MODIFY);

  npair = 0;
  double press = 0.0;
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

    ++npair;

    press = comp_pressure->compute_scalar() * domain->xprd * domain->yprd * domain->zprd;
    press += comp_ke->compute_scalar();
    press /= 3.0 * vol;

    timer->stamp(Timer::MODIFY);

  } while(fabs(press) > pairtol && ((pairmax > 0) ? npair < pairmax : true));

  return press;
}

/* ---------------------------------------------------------------------- */

void NufebRun::module_reactor()
{
  // reset to biological timestep

  update->dt = biodt;
  reset_dt();

  // call all fixes implementing reactor_nufeb()

  if (modify->n_reactor_nufeb) {
    timer->stamp();
    modify->reactor_nufeb();
    timer->stamp(Timer::MODIFY);
  }
}

/* ---------------------------------------------------------------------- */

double NufebRun::get_time()
{
  double rv = 0.0;

#ifdef _WIN32

  // from MSD docs.
  FILETIME ct,et,kt,ut;
  union { FILETIME ft; uint64_t ui; } cpu;
  if (GetProcessTimes(GetCurrentProcess(),&ct,&et,&kt,&ut)) {
    cpu.ft = ut;
    rv = cpu.ui * 0.0000001;
  }

#else /* ! _WIN32 */

  struct rusage ru;
  if (getrusage(RUSAGE_SELF, &ru) == 0) {
    rv = (double) ru.ru_utime.tv_sec;
    rv += (double) ru.ru_utime.tv_usec * 0.000001;
  }

#endif /* ! _WIN32 */

  return rv;
}
