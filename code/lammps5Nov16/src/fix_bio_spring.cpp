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

/* ----------------------------------------------------------------------

------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_bio_spring.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "error.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1.0e-10

enum{TETHER,COUPLE};

/* ---------------------------------------------------------------------- */

FixBioSpring::FixBioSpring(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  group2(NULL)
{
  if (narg < 9) error->all(FLERR,"Illegal fix bio/spring  command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  dynamic_group_allow = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  xc = yc = zc = 0;
  tflag = 0;

  var = new char*[4];
  ivar = new int[4];

  if (strcmp(arg[3],"tether") == 0) {
    if (narg != 13) error->all(FLERR,"Illegal fix bio/spring  command");
    styleflag = TETHER;
    k_spring = force->numeric(FLERR,arg[4]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[5],"NULL") == 0) xflag = 0;
    else xc = force->numeric(FLERR,arg[5]);
    if (strcmp(arg[6],"NULL") == 0) yflag = 0;
    else yc = force->numeric(FLERR,arg[6]);
    if (strcmp(arg[7],"NULL") == 0) zflag = 0;
    else zc = force->numeric(FLERR,arg[7]);
    r0 = force->numeric(FLERR,arg[8]);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix bio/spring  command");

    for (int i = 0; i < 4; i++) {
      int n = strlen(&arg[9+i][2]) + 1;
      var[i] = new char[n];
      strcpy(var[i], &arg[9+i][2]);
    }

  } else if (strcmp(arg[3],"couple") == 0) {
    if (narg != 10) error->all(FLERR,"Illegal fix bio/spring  command");
    styleflag = COUPLE;

    int n = strlen(arg[4]) + 1;
    group2 = new char[n];
    strcpy(group2,arg[4]);
    igroup2 = group->find(arg[4]);
    if (igroup2 == -1)
      error->all(FLERR,"Fix bio/spring  couple group ID does not exist");
    if (igroup2 == igroup)
      error->all(FLERR,"Two groups cannot be the same in fix bio/spring  couple");
    group2bit = group->bitmask[igroup2];

    k_spring = force->numeric(FLERR,arg[5]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[6],"NULL") == 0) xflag = 0;
    else xc = force->numeric(FLERR,arg[6]);
    if (strcmp(arg[7],"NULL") == 0) yflag = 0;
    else yc = force->numeric(FLERR,arg[7]);
    if (strcmp(arg[8],"NULL") == 0) zflag = 0;
    else zc = force->numeric(FLERR,arg[8]);
    r0 = force->numeric(FLERR,arg[9]);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix bio/spring  command");

  } else error->all(FLERR,"Illegal fix spring command");

  ftotal[0] = ftotal[1] = ftotal[2] = ftotal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixBioSpring::~FixBioSpring()
{
  delete [] group2;

  int i;
  for (i = 0; i < 4; i++) {
    delete[] var[i];
  }

  delete[] var;
  delete[] ivar;
}

/* ---------------------------------------------------------------------- */

int FixBioSpring::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::init()
{
  // recheck that group 2 has not been deleted

  if (group2) {
    igroup2 = group->find(group2);
    if (igroup2 == -1)
      error->all(FLERR,"Fix bio/spring  couple group ID does not exist");
    group2bit = group->bitmask[igroup2];
  }

  masstotal = group->mass(igroup);
  if (styleflag == COUPLE) masstotal2 = group->mass(igroup2);

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }

  if (styleflag == TETHER) {
    for (int n = 0; n < 4; n++) {
      ivar[n] = input->variable->find(var[n]);
      if (ivar[n] < 0)
        error->all(FLERR, "Variable name for fix bio/spring does not exist");
      if (!input->variable->equalstyle(ivar[n]))
        error->all(FLERR, "Variable for fix bio/spring  is invalid style");
    }

    cd = input->variable->compute_equal(ivar[0]);
    twitch = input->variable->compute_equal(ivar[1]);
    motility = input->variable->compute_equal(ivar[2]);
    detach = input->variable->compute_equal(ivar[3]);
  }
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::post_force(int vflag)
{
  if (styleflag == TETHER) spring_tether();
  else spring_couple();
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::spring_tether()
{
  double xcm[3];

  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);

  group->xcm(igroup,masstotal,xcm);

  double dx,dy,dz,fx,fy,fz,r,dr;

  // tethering point (x,y,z) is decided on the y coordinate of the center of mass of the group

  if (!xflag && !yflag && !zflag && !tflag && xcm[1]<cd) {
    xc = xcm[0];
    yc = 0;
    zc = xcm[2];
    // tethering point will not be changed under cd condition
    tflag = 1;
  }

  // update the tethering point along the x axis
  if (twitch < 1)
    xc = xc + motility*update->dt;

  dx = xcm[0] - xc;
  dy = xcm[1] - yc;
  dz = xcm[2] - zc;
//  if (!xflag) dx = 0.0;
//  if (!yflag) dy = 0.0;
//  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  r = MAX(r,SMALL);
  dr = r - r0;

  fx = k_spring*dx*dr/r;
  fy = k_spring*dy*dr/r;
  fz = k_spring*dz*dr/r;
  ftotal[0] = -fx;
  ftotal[1] = -fy;
  ftotal[2] = -fz;
  ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
  if (dr < 0.0) ftotal[3] = -ftotal[3];
  espring = 0.5*k_spring * dr*dr;

  if (masstotal > 0.0) {
    fx /= masstotal;
    fy /= masstotal;
    fz /= masstotal;
  }

  // if the distance between the centre of atoms and tethering point is greater than detach,
  // no spring force updated
  if (r > detach)
    return;

  // apply restoring force to atoms in group

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massone;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] -= fx*massone;
        f[i][1] -= fy*massone;
        f[i][2] -= fz*massone;
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] -= fx*massone;
        f[i][1] -= fy*massone;
        f[i][2] -= fz*massone;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::spring_couple()
{
  double xcm[3],xcm2[3];

  if (group->dynamic[igroup])
    masstotal = group->mass(igroup);

  if (group->dynamic[igroup2])
    masstotal2 = group->mass(igroup2);

  group->xcm(igroup,masstotal,xcm);
  group->xcm(igroup2,masstotal2,xcm2);

  // fx,fy,fz = components of k * (r-r0) / masstotal
  // fx2,fy2,fz2 = components of k * (r-r0) / masstotal2

  double dx,dy,dz,fx,fy,fz,fx2,fy2,fz2,r,dr;

  dx = xcm2[0] - xcm[0] - xc;
  dy = xcm2[1] - xcm[1] - yc;
  dz = xcm2[2] - xcm[2] - zc;
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  r = MAX(r,SMALL);
  dr = r - r0;

  fx = k_spring*dx*dr/r;
  fy = k_spring*dy*dr/r;
  fz = k_spring*dz*dr/r;
  ftotal[0] = fx;
  ftotal[1] = fy;
  ftotal[2] = fz;
  ftotal[3] = sqrt(fx*fx + fy*fy + fz*fz);
  if (dr < 0.0) ftotal[3] = -ftotal[3];
  espring = 0.5*k_spring * dr*dr;

  if (masstotal2 > 0.0) {
    fx2 = fx/masstotal2;
    fy2 = fy/masstotal2;
    fz2 = fz/masstotal2;
  } else fx2 = fy2 = fz2 = 0.0;

  if (masstotal > 0.0) {
    fx /= masstotal;
    fy /= masstotal;
    fz /= masstotal;
  } else fx = fy = fz = 0.0;

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double massone;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massone = rmass[i];
        f[i][0] += fx*massone;
        f[i][1] += fy*massone;
        f[i][2] += fz*massone;
      }
      if (mask[i] & group2bit) {
        massone = rmass[i];
        f[i][0] -= fx2*massone;
        f[i][1] -= fy2*massone;
        f[i][2] -= fz2*massone;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massone = mass[type[i]];
        f[i][0] += fx*massone;
        f[i][1] += fy*massone;
        f[i][2] += fz*massone;
      }
      if (mask[i] & group2bit) {
        massone = mass[type[i]];
        f[i][0] -= fx2*massone;
        f[i][1] -= fy2*massone;
        f[i][2] -= fz2*massone;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBioSpring::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of stretched spring
------------------------------------------------------------------------- */

double FixBioSpring::compute_scalar()
{
  return espring;
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixBioSpring::compute_vector(int n)
{
  return ftotal[n];
}
