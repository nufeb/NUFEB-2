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

#include "fix_adhesion_bacillus.h"

#include <math.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "error.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "group.h"
#include "pair.h"
#include "force.h"
#include "atom_vec_bacillus.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1e-30

/* ---------------------------------------------------------------------- */

FixAdhesionBacillus::FixAdhesionBacillus(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal fix nufeb/adhesion/eps command");

  ke = utils::numeric(FLERR,arg[3],true,lmp);
  cutoff = utils::numeric(FLERR,arg[4],true,lmp);
  if (cutoff < 0) error->all(FLERR, "Illegal value for cutoff parameter in fix nufeb/adhesion/bacillus");

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"Pair bacillus requires "
                        "atom style bacillus");
}

/* ---------------------------------------------------------------------- */

int FixAdhesionBacillus::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdhesionBacillus::init() {
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->size = 1;
}

/* ---------------------------------------------------------------------- */

void FixAdhesionBacillus::init_list(int id, NeighList *ptr) {
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixAdhesionBacillus::post_force(int vflag)
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixAdhesionBacillus::compute()
{
  int i,j,ii,jj,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq, leni, lenj, radi, radj;
  int *jlist;
  
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  int *bacillus = atom->bacillus;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *rmass = atom->rmass;

  AtomVecBacillus::Bonus *ibonus;
  AtomVecBacillus::Bonus *jbonus;

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      if (bacillus[i] >= 0) {
	int index = atom->bacillus[i];
	ibonus = &avec->bonus[index];
	leni = ibonus->length/2;
	radi = ibonus->diameter/2;
      }

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	j &= NEIGHMASK;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (bacillus[i] < 0 || bacillus[j] < 0) continue;

	int index = atom->bacillus[j];
	jbonus = &avec->bonus[index];
	lenj = jbonus->length/2;
	radj = jbonus->diameter/2;

	// no interaction
	double r = sqrt(rsq);
	if (r > radi+radj+leni+lenj+cutoff) continue;

	// rod-rod interaction
	rod_against_rod(i, j, x, v, f, torque, angmom, rmass, ibonus, jbonus);
      }
    }
  }
}

/* ----------------------------------------------------------------------*/
void FixAdhesionBacillus::rod_against_rod(int i, int j, double** x, double** v,
			           double** f, double** torque, double** angmom, double* rmass,
				   AtomVecBacillus::Bonus *&ibonus,
				   AtomVecBacillus::Bonus *&jbonus)
{
  double xi1[3],xi2[3],xpj1[3],xpj2[3];
  double r,t1,t2,h1[3],h2[3];
  double contact_dist, energy;

  avec->get_pole_coords(i, xi1, xi2);
  avec->get_pole_coords(j, xpj1, xpj2);

  contact_dist = (ibonus->diameter + jbonus->diameter)/2;

  int jflag = 1;

  distance_bt_rods(xpj1, xpj2, xi1, xi2, h2, h1, t2, t1, r);

  // include the vertices for interactions
  if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1 &&
      r < contact_dist + cutoff) {

    adhesion_force_and_torque(j, i, h2, h1, r, contact_dist, x, v, f, torque, angmom, rmass);
  }
}


/* ----------------------------------------------------------------------
  Compute adhesion forces and torques between two rods
------------------------------------------------------------------------- */

void FixAdhesionBacillus::adhesion_force_and_torque(int i, int j,
                 double* pi, double* pj, double r, double contact_dist,
		 double** x, double** v, double** f, double** torque,
		 double** angmom, double *rmass)
{
  double delx,dely,delz,R,fx,fy,fz;
  double del, rinv, ccel;
  double masssum;
  int newton_pair = force->newton_pair;

  delx = pi[0] - pj[0];
  dely = pi[1] - pj[1];
  delz = pi[2] - pj[2];
  R = r - contact_dist;
  masssum = rmass[i] + rmass[j];

  if (R > 0) {
    rinv = 1 / r;
    ccel = -masssum * ke * R;

    fx = delx * ccel * rinv;
    fy = dely * ccel * rinv;
    fz = delz * ccel * rinv;

    f[i][0] += fx;
    f[i][1] += fy;
    f[i][2] += fz;
    sum_torque(x[i], pi, fx, fy, fz, torque[i]);

    if (newton_pair) {
      f[j][0] -= fx;
      f[j][1] -= fy;
      f[j][2] -= fz;
      sum_torque(x[j], pj, -fx, -fy, -fz, torque[j]);
    }
  }
}

/* ----------------------------------------------------------------------
 compute the shortest distance between two rods (line segments)
------------------------------------------------------------------------- */

void FixAdhesionBacillus::distance_bt_rods(const double* x1,
		  const double* x2, const double* x3, const double* x4,
		  double* h1, double* h2, double& t1, double& t2, double& r)
{
    double ux = x2[0] - x1[0];
    double uy = x2[1] - x1[1];
    double uz = x2[2] - x1[2];

    double vx = x4[0] - x3[0];
    double vy = x4[1] - x3[1];
    double vz = x4[2] - x3[2];

    double wx = x1[0] - x3[0];
    double wy = x1[1] - x3[1];
    double wz = x1[2] - x3[2];

    double a = dot(ux, uy, uz, ux, uy, uz);         // always >= 0
    double b = dot(ux, uy, uz, vx, vy, vz);
    double c = dot(vx, vy, vz, vx, vy, vz);        // always >= 0
    double d = dot(ux, uy, uz, wx, wy, wz);
    double e = dot(vx, vy, vz, wx, wy, wz);
    double dd = a * c - b * b;
    double nt1, dt1 = dd;
    double nt2, dt2 = dd;

    // compute the line parameters of the two closest points
    if (dd < EPSILON) { // the lines are almost parallel
      nt1 = 0.0;
      dt1 = 1.0;         // to prevent possible division by 0.0 later
      nt2 = e;
      dt2 = c;
    } else {
      nt1 = (b*e - c*d);
      nt2 = (a*e - b*d);
      if (nt1 < 0.0) {
	nt1 = 0.0;
	nt2 = e;
	dt2 = c;
      } else if (nt1 > dt1) {
	nt1 = dt1;
	nt2 = e + b;
	dt2 = c;
      }
    }

    if (nt2 < 0.0) {
      nt2 = 0.0;

      if (-d < 0.0)
	nt1 = 0.0;
      else if (-d > a)
	nt1 = dt1;
      else {
	nt1 = -d;
	dt1 = a;
      }
    }
    else if (nt2 > dt2) {
      nt2 = dt2;
      // recompute sc for this edge
      if ((-d + b) < 0.0)
	nt1 = 0;
      else if ((-d + b) > a)
	nt1 = dt1;
      else {
	nt1 = (-d +  b);
	dt1 = a;
      }
    }

    // finally do the division to get t2 and t1
    t1 = (fabs(nt1) < EPSILON ? 0.0 : nt1 / dt1);
    t2 = (fabs(nt2) < EPSILON ? 0.0 : nt2 / dt2);

    double r1 = wx + t1*ux - t2*vx;
    double r2 = wy + t1*uy - t2*vy;
    double r3 = wz + t1*uz - t2*vz;

    r = norm(r1, r2, r3);   // return the closest distance

    h1[0] = x1[0] + (x2[0] - x1[0]) * t1;
    h1[1] = x1[1] + (x2[1] - x1[1]) * t1;
    h1[2] = x1[2] + (x2[2] - x1[2]) * t1;

    h2[0] = x3[0] + (x4[0] - x3[0]) * t2;
    h2[1] = x3[1] + (x4[1] - x3[1]) * t2;
    h2[2] = x3[2] + (x4[2] - x3[2]) * t2;
}

/* ----------------------------------------------------------------------
  Accumulate torque to rod from the force f=(fx,fy,fz) acting at point x
------------------------------------------------------------------------- */

void FixAdhesionBacillus::sum_torque(double* xm, double *x, double fx,
                                      double fy, double fz, double* torque)
{
  double rx = x[0] - xm[0];
  double ry = x[1] - xm[1];
  double rz = x[2] - xm[2];
  double tx = ry * fz - rz * fy;
  double ty = rz * fx - rx * fz;
  double tz = rx * fy - ry * fx;
  torque[0] += tx;
  torque[1] += ty;
  torque[2] += tz;
}
