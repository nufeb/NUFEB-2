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
#include "test_comm_grid.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "error.h"
#include "grid.h"
#include "comm_grid.h"
#include "grid_vec.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

TestCommGrid::TestCommGrid(LAMMPS *lmp, int narg, char **arg) :
  Integrate(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void TestCommGrid::init()
{
  Integrate::init();
}

/* ----------------------------------------------------------------------
   setup before run
------------------------------------------------------------------------- */

void TestCommGrid::setup(int flag)
{
  if (comm->me == 0 && screen) {
    fprintf(screen,"Setting up grid communication test run ...\n");
  }

  if (grid->nsubs < 1)
    error->all(FLERR, "Run style test/comm_grid requires at least one substrate");
  
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  atom->setup();
  domain->pbc();
  domain->reset_box();
  comm->setup();
  if (neighbor->style) neighbor->setup_bins();
  comm->exchange();
  comm->borders();
  domain->image_check();
  domain->box_too_small_check();
  neighbor->build(1);
  neighbor->ncalls = 0;
  grid->setup();
  comm_grid->setup();

  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
------------------------------------------------------------------------- */

void TestCommGrid::setup_minimal(int flag)
{
  update->setupflag = 1;

  // setup domain, communication and neighboring
  // acquire ghosts
  // build neighbor lists

  if (flag) {
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    comm->borders();
    domain->image_check();
    domain->box_too_small_check();
    neighbor->build(1);
    neighbor->ncalls = 0;
    grid->setup();
    comm_grid->setup();
  }

  update->setupflag = 0;
}

/* ----------------------------------------------------------------------
   run test
------------------------------------------------------------------------- */

void TestCommGrid::run(int)
{
  if (comm->me == 0 && screen) {
    fprintf(screen,"Running grid communication test run ...\n");
  }

  int *mask = grid->mask;
  double **conc = grid->conc;

  if (!conc)
    error->all(FLERR, "Run style comm_grid/test requires a grid style that uses conc attribute");
  
  for (int i = 0; i < grid->ncells; i++)
    conc[0][i] = 0;

  for (int z = grid->sublo[2] + 1; z < grid->subhi[2] - 1; z++) {
    for (int y = grid->sublo[1] + 1; y < grid->subhi[1] - 1; y++) {
      for (int x = grid->sublo[0] + 1; x < grid->subhi[0] - 1; x++) {
  	int i = (x - grid->sublo[0]) + (y - grid->sublo[1]) * grid->subbox[0] +
  	  (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
  	conc[0][i] = (x + 1) + (y + 1) * grid->extbox[0] +
  	  (z + 1) * grid->extbox[0] * grid->extbox[1];
      }
    }
  }

  if (comm->me == 0 && screen) {
    fprintf(screen,"Performing forward grid communication test ...");
  }

  const int nforward = 100;
  
  // call forward_comm a few times to stress it
  for (int i = 0; i < nforward; i++)
    comm_grid->forward_comm();

  // check for errors
  bool wrong = check();
  if (comm->me == 0 && screen) {
    if (wrong)
      fprintf(screen, " failed\n");
    else
      fprintf(screen, " passed\n");
  }

  if (comm->me == 0 && screen) {
    fprintf(screen,"Performing migration grid communication test ...");
  }

  // generate a test split modyfing the current domain decomposition
  double (*mysplit)[2] = comm->mysplit;
  double splitcpy[3][2];
  memcpy(splitcpy, mysplit, 6*sizeof(int));
  mysplit[0][0] = 0.0;
  mysplit[0][1] = 1.0;
  mysplit[1][0] = 0.0;
  mysplit[1][1] = 1.0;
  mysplit[2][0] = 0.0;
  mysplit[2][1] = 1.0;
  int splitaxis = 0;
  if (comm->nprocs > 1) {
    for (int i = 0; i <= comm->me; i++) {
      if (comm->me == comm->nprocs - 1 && i == comm->me) break;
      if (comm->me > i)
	mysplit[splitaxis][0] =
	  0.25 * (mysplit[splitaxis][1] - mysplit[splitaxis][0]);
      else mysplit[splitaxis][1] =
	     0.25 * (mysplit[splitaxis][1] - mysplit[splitaxis][0]);
      splitaxis = ++splitaxis % 3;
    }
  }
  int layout = comm->layout;
  comm->layout = Comm::LAYOUT_TILED;
  domain->set_local_box();
  comm_grid->migrate();
  comm_grid->setup();
  for (int i = 0; i < nforward; i++)
    comm_grid->forward_comm();
  
  wrong = check();
  if (!wrong && comm->me == 0 && screen) {
    if (wrong)
      fprintf(screen, " failed\n");
    else
      fprintf(screen, " passed\n");
  }

  comm->layout = layout;
  memcpy(mysplit, splitcpy, 6*sizeof(int));
  domain->set_local_box();
  grid->setup();
  comm_grid->setup();
}

bool TestCommGrid::check()
{
  bool result = false;
  const double small = 1e-12;
  int *mask = grid->mask;
  double **conc = grid->conc;
  for (int z = grid->sublo[2]; z < grid->subhi[2]; z++) {
    for (int y = grid->sublo[1]; y < grid->subhi[1]; y++) {
      for (int x = grid->sublo[0]; x < grid->subhi[0]; x++) {
	int i = (x - grid->sublo[0]) +
	  (y - grid->sublo[1]) * grid->subbox[0] +
	  (z - grid->sublo[2]) * grid->subbox[0] * grid->subbox[1];
	int j = (x + 1) + (y + 1) * grid->extbox[0] +
	  (z + 1) * grid->extbox[0] * grid->extbox[1];
	if (mask[i] & CORNER_MASK) {
	  continue;
	} else if ((mask[i] & X_NB_MASK) && domain->xperiodic) {
	  if (fabs(conc[0][i] - (j + grid->box[0])) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (-x periodic boundary): expected %d got %e\n",
		    comm->me, i, j + grid->box[0], conc[0][i]);
	  }
	} else if ((mask[i] & X_PB_MASK) && domain->xperiodic) {
	  if (fabs(conc[0][i] - (j - grid->box[0])) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (+x periodic boundary): expected %d got %e\n",
		    comm->me, i, j - grid->box[0], conc[0][i]);
	  }
	} else if ((mask[i] & Y_NB_MASK) && domain->yperiodic) {
	  if (fabs(conc[0][i] - (j + grid->extbox[0] * grid->box[1])) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (-y periodic boundary): expected %d got %e\n",
		    comm->me, i, j + grid->extbox[0] * grid->box[1], conc[0][i]);
	  }
	} else if ((mask[i] & Y_PB_MASK) && domain->yperiodic) {
	  if (fabs(conc[0][i] - (j - grid->extbox[0] * grid->box[1])) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (+y periodic boundary): expected %d got %e\n",
		    comm->me, i, j - grid->extbox[0] * grid->box[1], conc[0][i]);
	  }
	} else if ((mask[i] & Z_NB_MASK) && domain->zperiodic) {
	  if (fabs(conc[0][i] - (j + grid->extbox[0] * grid->extbox[1] * grid->box[2])) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (-z periodic boundary): expected %d got %e\n",
		    comm->me, i, j + grid->extbox[0] * grid->extbox[1] * grid->box[2], conc[0][i]);
	  }
	} else if ((mask[i] & Z_PB_MASK) && domain->zperiodic) {
	  if (fabs(conc[0][i] - (j - grid->extbox[0] * grid->extbox[1] * grid->box[2])) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (+z periodic boundary): expected %d got %e\n",
		    comm->me, i, j - grid->extbox[0] * grid->extbox[1] * grid->box[2], conc[0][i]);
	  }
	} else if (!(mask[i] & GHOST_MASK)) {
	  if (fabs(conc[0][i] - j) > small) {
	    result = true;
	    fprintf(screen, "[%d] Wrong value at cell %d (local domain): expected %d got %e\n",
		    comm->me, i, j, conc[0][i]);
	  }
	}
      }
    }
  }
  return result;
}
