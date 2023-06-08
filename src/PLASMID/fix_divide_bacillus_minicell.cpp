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

#include "fix_divide_bacillus_minicell.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "lmptype.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"
#include "compute.h"
#include "group.h"
#include "modify.h"
#include "domain.h"
#include "atom_masks.h"
#include "atom_vec_bacillus.h"
#include "random_park.h"

#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{SIZER,ADDER};
enum{MASS,LENGTH};

/* ---------------------------------------------------------------------- */

FixDivideBacillusMinicell::FixDivideBacillusMinicell(LAMMPS *lmp, int narg, char **arg) :
  FixDivide(lmp, narg, arg), birth_length(nullptr)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Fix nufeb/divide/bacillus/minicell requires "
      "atom style bacillus");

  restart_peratom = 1;

  if (narg < 8)
    error->all(FLERR, "Illegal fix nufeb/divide/bacillus/minicell command");

  imini = group->find(arg[3]);

  if (imini < 0)
    error->all(FLERR, "Can't find minicell group name");

  type = utils::inumeric(FLERR,arg[4],true,lmp);

  if (strcmp(arg[5], "sizer") == 0) {
    divflag = SIZER;
  } else if (strcmp(arg[5], "adder") == 0) {
    divflag = ADDER;
  } else {
    error->all(FLERR, "Illegal fix nufeb/division/bacillus/minicell command");
  }

  maxlength = utils::numeric(FLERR,arg[6],true,lmp);
  if (maxlength <= 0)
    error->all(FLERR, "Critical division length cannot be less or equal to 0");
  prob = utils::numeric(FLERR,arg[7],true,lmp);

  if (prob < 0 || prob > 1)
    error->all(FLERR, "Illegal minicell division probability value");
  seed = utils::numeric(FLERR,arg[8],true,lmp);

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  conserveflag = 0;
  int iarg = 9;
  var = 0.0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "normal") == 0) {
      var = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "conserve") == 0) {
      if (strcmp(arg[iarg+1], "mass") == 0) {
        conserveflag = MASS;
      } else if (strcmp(arg[iarg+1], "length") == 0) {
        conserveflag = LENGTH;
      } else {
        error->all(FLERR, "Illegal fix nufeb/division/bacillus/minicell command");
      }
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/division/bacillus/minicell command");
    }
  }

  maxradius = 0.0;

  // perform initial allocation of atom-based arrays
  // register with atom class

  peratom_flag = 1;
  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
}

FixDivideBacillusMinicell::~FixDivideBacillusMinicell()
{
  delete random;
  memory->destroy(birth_length);
  atom->delete_callback(id,Atom::GROW);
};

/* ---------------------------------------------------------------------- */

void FixDivideBacillusMinicell::init()
{
  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

      birth_length[i] = bonus->length;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixDivideBacillusMinicell::compute()
{  
  int nlocal = atom->nlocal;
  int mini_mask = group->bitmask[imini];
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      int ibonus = atom->bacillus[i];
      AtomVecBacillus::Bonus *bonus = &avec->bonus[ibonus];

      if (divflag == SIZER) divlength = maxlength;
      else if (divflag == ADDER) divlength = birth_length[i] + maxlength;
      // calculate a gaussian number with mean=maxlength, sd=var
      if (var) divlength += var*random->gaussian();

      if (bonus->length < divlength) continue;

      double imass, jmass;
      double ilen, xp1[3], xp2[3];
      int j;

      double vsphere = four_thirds_pi * atom->radius[i]*atom->radius[i]*atom->radius[i];
      double acircle = MY_PI*atom->radius[i]*atom->radius[i];
      double density = atom->rmass[i] / (vsphere + acircle * bonus->length);

      double oldx = atom->x[i][0];
      double oldy = atom->x[i][1];
      double oldz = atom->x[i][2];
      double old_len = bonus->length;

      xp1[0] = xp1[1] = xp1[2] = 0.0;
      xp2[0] = xp2[1] = xp2[2] = 0.0;

      avec->get_pole_coords(i, xp1, xp2);

      double prob_mini = random->uniform();
      // normal division
      if (prob_mini > prob) {

        // conserve mass
        if (conserveflag == MASS) {
          imass = atom->rmass[i]/2;
          ilen = (imass / density - vsphere) / acircle;
        // conserve length
        } else {
          ilen = old_len/2;
          imass = density * (vsphere + acircle*ilen);
        }

        jmass = imass;

            // update parent
        double dl = (0.5*ilen + atom->radius[i]) / (0.5*old_len);
        atom->x[i][0] += (xp1[0] - oldx) * dl;
        atom->x[i][1] += (xp1[1] - oldy) * dl;
        atom->x[i][2] += (xp1[2] - oldz) * dl;

        atom->rmass[i] = imass;
        birth_length[i] = ilen;

        bonus->pole1[0] *= ilen / old_len;
        bonus->pole1[1] *= ilen / old_len;
        bonus->pole1[2] *= ilen / old_len;
        bonus->pole2[0] *= ilen / old_len;
        bonus->pole2[1] *= ilen / old_len;
        bonus->pole2[2] *= ilen / old_len;
        bonus->length = ilen;

        // create child
        double *coord = new double[3];

        coord[0] = oldx + (xp2[0] - oldx) * dl;
        coord[1] = oldy + (xp2[1] - oldy) * dl;
        coord[2] = oldz + (xp2[2] - oldz) * dl;

        avec->create_atom(atom->type[i], coord);
        j = atom->nlocal - 1;
        atom->bacillus[j] = 0;
        atom->mask[j] = atom->mask[i];
        avec->set_bonus(j, bonus->pole1, bonus->diameter, bonus->quat, bonus->inertia);
        double d = sqrt(bonus->pole1[0]*bonus->pole1[0] +
			bonus->pole1[1]*bonus->pole1[1] +
			bonus->pole1[2]*bonus->pole1[2]);
        birth_length[j] = d*2;
        delete[] coord;
      } else {
        // abnormal division, generate one sphere (j) and one long rod (i)
        double prob_pole = random->uniform();

        double *coord = new double[3];
        double idl, jdl, pole1[3];

        jmass = vsphere * density;
            imass = atom->rmass[i] - jmass;

            ilen = (imass / density - vsphere) / acircle;
            idl = ilen / old_len;
        jdl = (ilen + 2*atom->radius[i]) / old_len;

        atom->rmass[i] = imass;
            birth_length[i] = ilen;
        bonus->length = ilen;
        bonus->pole1[0] *= ilen / old_len;
        bonus->pole1[1] *= ilen / old_len;
        bonus->pole1[2] *= ilen / old_len;
        bonus->pole2[0] *= ilen / old_len;
        bonus->pole2[1] *= ilen / old_len;
        bonus->pole2[2] *= ilen / old_len;

        if (prob_pole > 0.5) {
          // update parent
          atom->x[i][0] = xp1[0] + (oldx - xp1[0]) * idl;
          atom->x[i][1] = xp1[1] + (oldy - xp1[1]) * idl;
          atom->x[i][2] = xp1[2] + (oldz - xp1[2]) * idl;

          // create child
          coord[0] = oldx + (xp2[0] - oldx) * jdl;
          coord[1] = oldy + (xp2[1] - oldy) * jdl;
          coord[2] = oldz + (xp2[2] - oldz) * jdl;
            } else {
          // update parent
          atom->x[i][0] = xp2[0] + (oldx - xp2[0]) * idl;
          atom->x[i][1] = xp2[1] + (oldy - xp2[1]) * idl;
          atom->x[i][2] = xp2[2] + (oldz - xp2[2]) * idl;

          // create child
          coord[0] = oldx + (xp1[0] - oldx) * jdl;
          coord[1] = oldy + (xp1[1] - oldy) * jdl;
          coord[2] = oldz + (xp1[2] - oldz) * jdl;
        }

        avec->create_atom(type, coord);
        j = atom->nlocal - 1;
        atom->bacillus[j] = 0;
        atom->mask[j] = 1 | mini_mask;
        pole1[0] = pole1[1] = pole1[2] = 0;
        avec->set_bonus(j, pole1, bonus->diameter, bonus->quat, bonus->inertia);
        birth_length[j] = 0.0;
        delete[] coord;
      }
      // set daughter j attributes
      atom->tag[j] = 0;
      atom->v[j][0] = atom->v[i][0];
      atom->v[j][1] = atom->v[i][1];
      atom->v[j][2] = atom->v[i][2];
      atom->f[j][0] = atom->f[i][0];
      atom->f[j][1] = atom->f[i][1];
      atom->f[j][2] = atom->f[i][2];
      atom->torque[j][0] = atom->torque[i][0];
      atom->torque[j][1] = atom->torque[i][1];
      atom->torque[j][2] = atom->torque[i][2];
      atom->angmom[j][0] = atom->angmom[i][0];
      atom->angmom[j][1] = atom->angmom[i][1];
      atom->angmom[j][2] = atom->angmom[i][2];
      atom->rmass[j] = jmass;
      atom->biomass[j] = atom->biomass[i];
      atom->radius[j] = atom->radius[i];

      modify->create_attribute(j);

      for (int m = 0; m < modify->nfix; m++) {
        modify->fix[m]->update_arrays(i, j);
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

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

void *FixDivideBacillusMinicell::extract(const char *str, int &itype)
{
  if (strcmp(str,"radius") == 0) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->bacillus[i] >= 0) {
        double radius = atom->radius[i];
        if (radius > maxradius) maxradius = radius;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &maxradius, 1, MPI_DOUBLE, MPI_MAX, world);

    return &maxradius;
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixDivideBacillusMinicell::grow_arrays(int nmax)
{
  memory->grow(birth_length,nmax,"fix nufeb/divide/bacillus/minicell:birth_length");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixDivideBacillusMinicell::copy_arrays(int i, int j, int /*delflag*/)
{
  birth_length[j] = birth_length[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixDivideBacillusMinicell::pack_exchange(int i, double *buf)
{
  buf[0] = birth_length[i];

  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixDivideBacillusMinicell::unpack_exchange(int nlocal, double *buf)
{
  birth_length[nlocal] = buf[0];
  return 1;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */
double FixDivideBacillusMinicell::memory_usage()
{
  double bytes = atom->nmax * sizeof(double);
  return bytes;
}

