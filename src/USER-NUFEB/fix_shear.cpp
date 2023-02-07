/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_shear.h"

#include <string.h>
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixShear::FixShear(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix nufeb/shear command");

  xflag = yflag = -1;
  layer = domain->prd[2];
  shear_rate = utils::numeric(FLERR,arg[3],false,lmp);
  viscosity = utils::numeric(FLERR,arg[4],false,lmp);


  if (strcmp(arg[5], "+x") == 0) {
    xflag = 1;
  } else if (strcmp(arg[5], "-x") == 0) {
    xflag = 0;
  } else if (strcmp(arg[5], "+y") == 0) {
    yflag = 1;
  } else if (strcmp(arg[5], "-y") == 0) {
    yflag = 0;
  } else {
    error->all(FLERR, "Illegal fix nufeb/shear command");
  }

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"layer") == 0) {
      double input;
      input = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if(input > layer || input < 0)
	 error->all(FLERR,"Illegal fix nufeb/shear command");
      layer = input;
      iarg += 2;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixShear::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixShear::post_force(int /*vflag*/)
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixShear::compute()
{
  double **f = atom->f;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  double *radius = atom->radius;

  int nlocal = atom->nlocal;
  double factor;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double diameter = 2 * radius[i];
      double my_3pi_vd = MY_3PI * viscosity * diameter;

      (x[i][2] / layer > 1) ? factor = 1 : factor = x[i][2] / layer;
      if (xflag == 1) {
	f[i][0] += my_3pi_vd * (shear_rate * factor - v[i][0]);
	f[i][1] += my_3pi_vd * (0.0 - v[i][1]);
	f[i][2] += my_3pi_vd * (0.0 - v[i][2]);
      } else if (xflag == 0){
	f[i][0] -= my_3pi_vd * (shear_rate * factor - v[i][0]);
	f[i][1] += my_3pi_vd * (0.0 - v[i][1]);
	f[i][2] += my_3pi_vd * (0.0 - v[i][2]);
      } else if (yflag == 1){
	f[i][0] += my_3pi_vd * (0.0 - v[i][0]);
	f[i][1] += my_3pi_vd * (shear_rate * factor - v[i][1]);
	f[i][2] += my_3pi_vd * (0.0 - v[i][2]);
      } else if (yflag == 0){
	f[i][0] += my_3pi_vd * (0.0 - v[i][0]);
	f[i][1] -= my_3pi_vd * (shear_rate * factor - v[i][1]);
	f[i][2] += my_3pi_vd * (0.0 - v[i][2]);
      }
    }
  }
}
