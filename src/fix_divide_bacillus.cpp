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

#include "fix_divide_bacillus.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_bacillus.h"
#include "error.h"
#include "force.h"
#include "lmptype.h"
#include "math_const.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define DELTA 1.005

/* ---------------------------------------------------------------------- */

FixDivideBacillus::FixDivideBacillus(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Fix nufeb/monod/ecoli/wild requires "
      "atom style bacillus");

  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/divide/bacillus command");
  
  maxlength = force->numeric(FLERR, arg[3]);
}

/* ---------------------------------------------------------------------- */

void FixDivideBacillus::compute()
{  
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

      if (bonus->length >= maxlength) {
	double *pole1 = bonus->pole1;
	double *pole2 = bonus->pole2;

	double parentx = atom->x[i][0];
	double parenty = atom->x[i][1];
	double parentz = atom->x[i][2];
	double parentpx = pole2[0];
	double parentpy = pole2[1];
	double parentyz = pole2[2];
	double shiftx = 2*(pole1[0] - parentx) / bonus->length * bonus->diameter;
	double shifty = 2*(pole1[1] - parenty) / bonus->length * bonus->diameter;
	double shiftz = 2*(pole1[2] - parentz) / bonus->length * bonus->diameter;

        double newx = (parentx + pole1[0])/2 + shiftx;
        double newy = (parenty + pole1[1])/2 + shifty;
        double newz = (parentz + pole1[2])/2 + shiftz;

        double new_rmass = atom->rmass[i]/2;
        double new_biomass = atom->biomass[i]/2;

        // update parent
        atom->rmass[i] = new_rmass;
        atom->biomass[i] = new_biomass;
        bonus->length /= 2;

        pole1[0] += shiftx;
        pole1[1] += shifty;
        pole1[2] += shiftz;
        pole2[0] = atom->x[i][0] + shiftx;
        pole2[1] = atom->x[i][1] + shifty;
        pole2[2] = atom->x[i][2] + shiftz;

        atom->x[i][0] = newx;
        atom->x[i][1] = newy;
        atom->x[i][2] = newz;

        // create child
        double *coord = new double[3];
        double *p1 = new double[3];
        double *p2 = new double[3];

        shiftx = (parentx - pole2[0]) / bonus->length * bonus->diameter / 2;
	shifty = (parenty - pole2[1]) / bonus->length * bonus->diameter / 2;
	shiftz = (parentz - pole2[2]) / bonus->length * bonus->diameter / 2;

        newx = (parentx + parentpx)/2 - shiftx;
        newy = (parenty + parentpy)/2 - shifty;
        newz = (parentz + parentyz)/2 - shiftz;

        coord[0] = newx;
        coord[1] = newy;
        coord[2] = newz;

        p1[0] = parentx - shiftx;
        p1[1] = parenty - shifty;
        p1[2] = parentz - shiftz;
        p2[0] = 2*newx - p1[0];
        p2[1] = 2*newy - p1[1];
        p2[2] = 2*newz - p1[2];

        avec->create_atom(atom->type[i], coord);
        int n = atom->nlocal - 1;
        atom->bacillus[n] = 0;
        avec->set_bonus(n, p1, p2, bonus);

        atom->tag[n] = 0;
        atom->mask[n] = atom->mask[i];
        atom->v[n][0] = atom->v[i][0];
        atom->v[n][1] = atom->v[i][1];
        atom->v[n][2] = atom->v[i][2];
	atom->f[n][0] = atom->f[i][0];
	atom->f[n][1] = atom->f[i][1];
	atom->f[n][2] = atom->f[i][2];
	atom->torque[n][0] = atom->torque[i][0];
	atom->torque[n][1] = atom->torque[i][1];
	atom->torque[n][2] = atom->torque[i][2];
	atom->angmom[n][0] = atom->angmom[i][0];
	atom->angmom[n][1] = atom->angmom[i][1];
	atom->angmom[n][2] = atom->angmom[i][2];
        atom->rmass[n] = new_rmass;
        atom->biomass[n] = new_biomass;
        atom->radius[n] = atom->radius[i];

        modify->create_attribute(n);

        delete[] coord;
        delete[] p1;
        delete[] p2;
      }
    }
  }
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
}
