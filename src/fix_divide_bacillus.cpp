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
#include "random_park.h"

#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixDivideBacillus::FixDivideBacillus(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Fix nufeb/monod/ecoli/wild requires "
      "atom style bacillus");

  if (narg < 5)
    error->all(FLERR, "Illegal fix nufeb/divide/bacillus command");
  
  maxlength = force->numeric(FLERR, arg[3]);
  seed = force->inumeric(FLERR, arg[4]);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
}

FixDivideBacillus::~FixDivideBacillus()
{
  delete random;
};

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
	double phix = random->uniform() * 2e-8;
	double phiy = random->uniform() * 2e-8;
	double phiz = random->uniform() * 2e-8;

	double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
	double acircle = MY_PI*atom->radius[i]*atom->radius[i];
	double density = atom->rmass[i] / (vsphere + acircle * bonus->length);

        double new_rmass = atom->rmass[i]/2;
        double new_biomass = atom->biomass[i]/2;

        double half_l = bonus->length / 2;
        // conserve mass
        double new_l = (new_rmass / density - vsphere) / acircle;
        double new_half_l = new_l / 2;
	double xp1[3];
	double xp2[3];

	xp1[0] = xp1[1] = xp1[2] = 0.0;
	xp2[0] = xp2[1] = xp2[2] = 0.0;

	avec->get_pole_coords(i, xp1, xp2, 0);

	double parentx = atom->x[i][0];
	double parenty = atom->x[i][1];
	double parentz = atom->x[i][2];
	double parent_length = bonus->length;

        // update parent
	double dl = (new_half_l + atom->radius[i]) / half_l;
	atom->x[i][0] += (xp1[0] - parentx) * dl;
	atom->x[i][1] += (xp1[1] - parenty) * dl;
	atom->x[i][2] += (xp1[2] - parentz) * dl;

        atom->rmass[i] = new_rmass;
        atom->biomass[i] = new_biomass;

        bonus->pole1[0] *= new_half_l/half_l;
        bonus->pole1[1] *= new_half_l/half_l;
        bonus->pole1[2] *= new_half_l/half_l;
        bonus->pole2[0] *= new_half_l/half_l;
        bonus->pole2[1] *= new_half_l/half_l;
        bonus->pole2[2] *= new_half_l/half_l;
        bonus->length = new_l;

//        printf("n: %i \n", i);
//	printf("c1=%e, c2=%e, c3=%e \n", atom->x[i][0],atom->x[i][1],atom->x[i][2]);
//	printf("p1x=%e, p1y=%e, p1z=%e \n", bonus->pole1[0],bonus->pole1[1],bonus->pole1[2]);
//	printf("p2x=%e, p2y=%e, p2z=%e \n", bonus->pole2[0],bonus->pole2[1],bonus->pole2[2]);
//	printf("\n");

        // create child
        double *coord = new double[3];

	coord[0] = parentx + (xp2[0] - parentx) * dl;
	coord[1] = parenty + (xp2[1] - parenty) * dl;
	coord[2] = parentz + (xp2[2] - parentz) * dl;

        avec->create_atom(atom->type[i], coord);
        int n = atom->nlocal - 1;
        atom->bacillus[n] = 0;

        avec->set_bonus(n, bonus->pole1, bonus->diameter, bonus->quat, bonus->inertia);
        ibonus = atom->bacillus[n];
        bonus = &avec->bonus[ibonus];
//        printf("n: %i \n", n);
//	printf("c1=%e, c2=%e, c3=%e \n", coord[0],coord[1],coord[2]);
//	printf("p1x=%e, p1y=%e, p1z=%e \n", bonus->pole1[0],bonus->pole1[1],bonus->pole1[2]);
//	printf("p2x=%e, p2y=%e, p2z=%e \n", bonus->pole2[0],bonus->pole2[1],bonus->pole2[2]);
//	printf("\n");
        double density2 = new_rmass / (vsphere + acircle * avec->bonus[atom->bacillus[n]].length);

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
