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
#include "intersect_list.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

IntersectList::IntersectList(LAMMPS *lmp, int m) : Pointers(lmp)
{
  n = 0;
  nmax = m;
  procs = memory->create(procs, m, "intersect_list:procs");
  boxlo = memory->create(boxlo, 3*m, "intersect_list:boxlo");
  boxhi = memory->create(boxhi, 3*m, "intersect_list:boxhi");
  count = memory->create(count, m, "intersect_list:count");
}

IntersectList::~IntersectList()
{
  memory->destroy(procs);
  memory->destroy(boxlo);
  memory->destroy(boxhi);
  memory->destroy(count);
}

void IntersectList::add(int proc, int *lo, int *hi, int c)
{
  if (nmax < n + 1) {
    nmax *= 2;
    procs = memory->grow(procs, nmax, "intersect_list:procs");
    boxlo = memory->grow(boxlo, 3*nmax, "intersect_list:boxlo");
    boxhi = memory->grow(boxhi, 3*nmax, "intersect_list:boxhi");
    count = memory->grow(count, nmax, "intersect_list:count");
  }
  procs[n] = proc;
  memcpy(&boxlo[3*n], lo, 3*sizeof(int));
  memcpy(&boxhi[3*n], hi, 3*sizeof(int));
  count[n] = c;
  ++n;
}
