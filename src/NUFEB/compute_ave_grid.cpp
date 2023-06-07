/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "compute_ave_grid.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "error.h"
#include "pointers.h"
#include "update.h"
#include "memory.h"
#include "group.h"

#include "modify.h"
#include "grid.h"
#include "grid_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeAveGrid::ComputeAveGrid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute nufeb/ave_grid command");

  vector_flag = 1;
  extvector = 0;

  gro_flag = rea_flag = den_flag = con_flag = 0;
  ph_flag = act_flag = 0;

  if (strcmp(arg[3], "gro") == 0) {
    gro_flag = 1;
    size_vector = 2;
  } else if (strcmp(arg[3], "rea") == 0) {
    rea_flag = 1;
    size_vector = grid->nsubs;
  } else if (strcmp(arg[3], "con") == 0) {
    con_flag = 1;
    size_vector = grid->nsubs;
  } else if (strcmp(arg[3], "den") == 0) {
    den_flag = 1;
    size_vector = 1;
  } else if (strcmp(arg[3], "ph") == 0) {
    ph_flag = 1;
    size_vector = 1;
  } else if (strcmp(arg[3], "act") == 0) {
    act_flag = 1;
    size_vector = grid->nsubs;
  } else
    error->all(FLERR,"Unknown keyword in compute nufeb/ave_grid");

  memory->create(vector,size_vector,"compute:vector");
}

/* ---------------------------------------------------------------------- */

ComputeAveGrid::~ComputeAveGrid()
{
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputeAveGrid::compute_vector()
{
  invoked_vector = update->ntimestep;

  if (act_flag) {
    if (grid->act == nullptr)
      error->all(FLERR, "Illegal compute nufeb/ave_grid command: act");
  } else if (ph_flag) {
    if (grid->ph == nullptr)
      error->all(FLERR,"Illegal compute nufeb/ave_grid command: ph");
  }

  for (int i = 0; i < size_vector; i++) {
    double sum = 0;
    int sum_cells = 0;
    for (int j = 0; j < grid->ncells; j++) {
      if (con_flag) {
        if (grid->mask[j] & GRID_MASK) {
          sum += grid->conc[i][j];
          sum_cells++;
        }
      } else if (rea_flag) {
        if (grid->mask[j] & GRID_MASK) {
          sum += grid->reac[i][j];
          sum_cells++;
        }
      } else if (ph_flag) {
        if (grid->mask[j] & GRID_MASK) {
          sum += grid->ph[j];
          sum_cells++;
        }
      } else if (act_flag) {
        if (grid->mask[j] & GRID_MASK) {
          sum += grid->act[i][j];
          sum_cells++;
        }
      } else if (den_flag) {
        if (grid->dens[igroup][j] > 0) {
          sum += grid->dens[igroup][j];
          sum_cells++;
        }
      } else if (gro_flag) {
        if (grid->growth[igroup][j][i] > 0) {
          sum += grid->growth[igroup][j][i];
          sum_cells++;
        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(MPI_IN_PLACE, &sum_cells, 1, MPI_INT, MPI_SUM, world);

    sum /= sum_cells;
    vector[i] = sum;
  }
}
