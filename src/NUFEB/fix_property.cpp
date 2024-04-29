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

#include "fix_property.h"

#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixProperty::FixProperty(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  aprop(nullptr), vprop(nullptr)
{
  restart_peratom = 1;
  // this fix produces either a per-atom vector or array
  peratom_flag = 1;
  // perform initial allocation of atom-based array
  // register with Atom class

  atom->add_callback(Atom::GROW);
  atom->add_callback(Atom::RESTART);
}

/* ---------------------------------------------------------------------- */

FixProperty::~FixProperty()
{
  if (copymode) return;

  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,Atom::GROW);
  atom->delete_callback(id,Atom::RESTART);

  memory->destroy(aprop);
  memory->destroy(vprop);
}

/* ---------------------------------------------------------------------- */

int FixProperty::modify_param(int narg, char **arg)
{
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "nevery") == 0) {
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nevery <= 0) error->all(FLERR,"Illegal fix_modify command");
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix_modify command");
    }
  }
  return iarg;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixProperty::memory_usage()
{
  double bytes;
  if (size_peratom_cols > 1)
    bytes = atom->nmax*size_peratom_cols*sizeof(double);
  else
    bytes = atom->nmax*sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixProperty::grow_arrays(int nmax)
{
  if (size_peratom_cols > 1) {
    memory->grow(aprop,nmax,size_peratom_cols,"fix_property:array");
    array_atom = aprop;
  } else {
    memory->grow(vprop,nmax,"fix_property:vector");
    vector_atom = vprop;
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixProperty::copy_arrays(int i, int j, int /*delflag*/)
{
  if (size_peratom_cols > 1) {
    for (int m = 0; m < size_peratom_cols; m++)
      aprop[j][m] = aprop[i][m];
  } else {
    vprop[j] = vprop[i];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixProperty::pack_exchange(int i, double *buf)
{
  if (size_peratom_cols > 1) {
    for (int m = 0; m < size_peratom_cols; m++)
      buf[m] = aprop[i][m];
  } else {
    buf[0] = vprop[i];
  }
  return size_peratom_cols + 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixProperty::unpack_exchange(int nlocal, double *buf)
{
  if (size_peratom_cols > 1) {
    for (int m = 0; m < size_peratom_cols; m++)
      aprop[nlocal][m] = buf[m];
  } else {
    vprop[nlocal] = buf[0];
  }
  return size_peratom_cols + 1;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixProperty::pack_restart(int i, double *buf)
{
  buf[0] = size_peratom_cols + 1;
  if (size_peratom_cols > 1) {
    for (int m = 0; m < size_peratom_cols; m++)
      buf[m+1] = aprop[i][m];
  } else {
    buf[1] = vprop[i];
  }
  return size_peratom_cols + 1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixProperty::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;
  if (size_peratom_cols > 1) {
    for (int i = 0; i < size_peratom_cols; i++)
      aprop[nlocal][i] = extra[nlocal][m++];
  } else {
    vprop[nlocal] = extra[nlocal][m];
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixProperty::maxsize_restart()
{
  return size_peratom_cols + 1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixProperty::size_restart(int /*nlocal*/)
{
  return size_peratom_cols + 1;
}

