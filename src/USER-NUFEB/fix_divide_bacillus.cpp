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
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

      if (bonus->length >= maxlength) {
	double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
	double acircle = MY_PI*atom->radius[i]*atom->radius[i];
	double density = atom->rmass[i] / (vsphere + acircle * bonus->length);

        double new_rmass = atom->rmass[i]/2;
        double new_biomass = atom->biomass[i]/2;

        double half_length = bonus->length / 2;
        // conserve mass
        double new_length = (new_rmass - density * vsphere) / (density * acircle);
        // conserve volume
        //double new_length = half_length - atom->radius[i];

	double *pole1 = bonus->pole1;
	double *pole2 = bonus->pole2;

	double parentx = atom->x[i][0];
	double parenty = atom->x[i][1];
	double parentz = atom->x[i][2];
	double parentpx = pole2[0];
	double parentpy = pole2[1];
	double parentpz = pole2[2];
	double parent_length = bonus->length;

        // update parent
	double dl = atom->radius[i] / half_length;
	pole2[0] = parentx - (parentx - pole1[0]) * dl;
	pole2[1] = parenty - (parenty - pole1[1]) * dl;
	pole2[2] = parentz - (parentz - pole1[2]) * dl;

	dl = (new_length - half_length) / half_length;
        pole1[0] -= (pole1[0] - parentx) * dl;
        pole1[1] -= (pole1[1] - parenty) * dl;
        pole1[2] -= (pole1[2] - parentz) * dl;

        atom->x[i][0] = (pole1[0] + pole2[0])/2;
        atom->x[i][1] = (pole1[1] + pole2[1])/2;;
        atom->x[i][2] = (pole1[2] + pole2[2])/2;;

        atom->rmass[i] = new_rmass;
        atom->biomass[i] = new_biomass;
        bonus->length = new_length;

        // create child
        double *coord = new double[3];
        double *p1 = new double[3];
        double *p2 = new double[3];

	dl = atom->radius[i] / half_length;
	p1[0] = parentx - (parentx - parentpx) * dl;
	p1[1] = parenty - (parenty - parentpy) * dl;
	p1[2] = parentz - (parentz - parentpz) * dl;

	dl = (new_length - half_length) / half_length;
        p2[0] = parentpx - (parentpx - parentx) * dl;
        p2[1] = parentpy - (parentpy - parenty) * dl;
        p2[2] = parentpz - (parentpz - parentz) * dl;

        // conserve volume
//        p2[0] = parentpx;
//        p2[1] = parentpy;
//        p2[2] = parentpz;

        coord[0] = (p1[0] + p2[0])/2;;
        coord[1] = (p1[1] + p2[1])/2;;
        coord[2] = (p1[2] + p2[2])/2;;

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
