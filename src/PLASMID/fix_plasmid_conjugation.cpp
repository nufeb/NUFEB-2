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

#include <string.h>
#include <math.h>
#include "atom.h"
#include "error.h"
#include "modify.h"
#include "group.h"
#include "update.h"
#include "memory.h"
#include "atom_masks.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "force.h"

#include "fix_plasmid_conjugation.h"
#include "fix_property_plasmid.h"
#include "atom_vec_bacillus.h"

#include "fix.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1e-30
#define PLM_TRANSFER_DIST 5e-7
#define CUTOFF 1e-7

/* ---------------------------------------------------------------------- */

FixPlasmidConjugation::FixPlasmidConjugation(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR, "Illegal nufeb/plasmid/conjugate command");

  fix_plm = nullptr;
  tflag = 0;

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for fix plasmid/conjugation");
  fix_plm = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  irecip = group->find(arg[3]);
  irecipbit = group->bitmask[irecip];
  itrans = group->find(arg[4]);
  itransbit = group->bitmask[itrans];

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"fix nufeb/plasmid/conjugate requires "
                        "atom style bacillus");
  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "type") == 0) {
      ttrans = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      tflag = 1;
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix nufeb/plasmid/conjugate command");
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixPlasmidConjugation::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPlasmidConjugation::biology_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixPlasmidConjugation::compute()
{
  int i,j,ii,jj,jnum;
  double xtmp,ytmp,ztmp;
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

  NeighList *list = force->pair->list;
  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(atom->mask[i] & groupbit)) continue;
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
      if (!(atom->mask[j] & irecipbit)) continue;

      double r, contact_dist;
      double xi1[3],xi2[3],xpj1[3],xpj2[3];
      double t1,t2,h1[3],h2[3];
      double delx,dely,delz;
      double R;

      double *vprop = fix_plm->vprop;
      double *plm_x = fix_plm->plm_x[i];
      double px[3], dmin;
      int plm;

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
      r = sqrt(rsq);
      if (r > radi+radj+leni+lenj) continue;

      avec->get_pole_coords(i, xi1, xi2);
      avec->get_pole_coords(j, xpj1, xpj2);

      contact_dist = (ibonus->diameter + jbonus->diameter)/2;
      distance_bt_rods(xpj1, xpj2, xi1, xi2, h2, h1, t2, t1, r);
      // rod-rod interaction
      R = r - contact_dist;
      dmin = 1e10;

      for (int m = 0; m < static_cast<int>(vprop[i]); m++) {
	fix_plm->get_plasmid_coords(i, m, px);
	delx = px[0] - h1[0];
	dely = px[1] - h1[1];
	delz = px[2] - h1[2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq < PLM_TRANSFER_DIST && rsq < dmin) {
	  dmin = rsq;
	  plm = m;
	}
      }

      if (R < CUTOFF) conjugate(j, h2);
    }
  }
}

/* ----------------------------------------------------------------------
   transfer plasmid m from bacillus i to bacillus j
------------------------------------------------------------------------- */
void FixPlasmidConjugation::conjugate(int j, double* h2) {
  double* vprop = fix_plm->vprop;
  double** nproteins = fix_plm->nproteins;

  int n = static_cast<int>(fix_plm->vprop[j]);
  if (fix_plm->rep_flag) nproteins[j][n] = 0;
  fix_plm->set_plm_x(j, n, h2, atom->x[j]);
  vprop[j]++;

  atom->mask[j] = itrans;
  if (tflag)
   atom->type[j] = ttrans;
}

/* ----------------------------------------------------------------------
   get quaternion from two vectors
------------------------------------------------------------------------- */
void FixPlasmidConjugation::get_quat(double *vec1, double *vec2, double *q){
  double ans1[3], ans2[3], cross[3];
  ans1[0] = ans1[1] = ans1[2] = 0.0;
  ans2[0] = ans2[1] = ans2[2] = 0.0;

  MathExtra::normalize3(vec1,ans2);
  MathExtra::normalize3(vec2,ans1);
  MathExtra::cross3(ans1, ans2, cross);

  double d = MathExtra::dot3(ans1, ans2);
  double s = sqrt((1+d)*2);
  double invs = 1 / s;

  q[0] = s*0.5;
  q[1] = cross[0]*invs;
  q[2] = cross[1]*invs;
  q[3] = cross[2]*invs;

  MathExtra::qnormalize(q);
}

/* ----------------------------------------------------------------------
 compute the shortest distance between two rods (line segments)
------------------------------------------------------------------------- */

void FixPlasmidConjugation::distance_bt_rods(const double* x1,
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
