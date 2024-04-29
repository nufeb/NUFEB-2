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

#include <cstdlib>
#include <cstring>
#include <math.h>
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "grid.h"
#include "atom_vec_bacillus.h"
#include "fix_divide_bacillus_minicell.h"
#include "update.h"
#include "math_const.h"
#include "math_extra.h"
#include "random_park.h"
#include <random>

#include "fix_property_plasmid.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyPlasmid::FixPropertyPlasmid(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg),
  plm_x(nullptr), fila(nullptr), pre_x(nullptr), nproteins(nullptr),
  tfila(nullptr), nfilas(nullptr)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"fix nufeb/property/plasmid requires "
      "atom style bacillus");

  if (narg < 3) error->all(FLERR,"Illegal fix nufeb/property/plasmid command");

  fila_max = 0;
  plm_max = 5;
  plm_init = 0;
  plm_dia = 5e-7;
  rep_flag = par_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "max") == 0) {
      plm_max = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      if (plm_max <= 0)
	error->all(FLERR,"Illegal fix nufeb/property/plasmid command: max");
      iarg += 2;
    } else if (strcmp(arg[iarg], "init") == 0) {
      plm_init = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      if (plm_init < 0)
	error->all(FLERR,"Illegal fix nufeb/property/plasmid command: init");
      iarg += 2;
    } else if (strcmp(arg[iarg], "dia") == 0) {
      plm_dia = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      if (plm_dia <= 0)
	error->all(FLERR,"Illegal fix nufeb/property/plasmid command: plm_dia");
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/property/plasmid command");
    }
  }

  if (plm_init > plm_max) error->all(FLERR,"Illegal fix nufeb/property/plasmid command: "
      "initial plasmid cannot be more than maximum plasmid number");

  size_peratom_cols = 1;
  fila_max = plm_max * (plm_max - 1) / 2;

  grow_arrays(atom->nmax);

  std::random_device rd;
  std::mt19937 eng(rd());
  for (int i = 0; i < atom->nlocal; i++) {
    vprop[i] = 0;
    if (!(atom->mask[i] & groupbit)) continue;
    AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
    vprop[i] = plm_init;
    // initialise plasmid position
    for (int j = 0; j < static_cast<int>(vprop[i]); j++) {
      double limit[6];

      int jx0 = j*3;
      int jx1 = j*3+1;
      int jx2 = j*3+2;

      // assign x
      limit[1] = (bouns->length + bouns->diameter - plm_dia) / 2;
      limit[0] = -limit[1];

      std::uniform_real_distribution<double> x(limit[0], limit[1]);
      plm_x[i][jx0] = x(eng);

      get_cell_boundary(limit, i, j);

      // assign y, z
      std::uniform_real_distribution<double> y(limit[2], limit[3]);
      std::uniform_real_distribution<double> z(limit[4], limit[5]);
      plm_x[i][jx1] = y(eng);
      plm_x[i][jx2] = z(eng);
    }
  }
}

/* ---------------------------------------------------------------------- */
FixPropertyPlasmid:: ~FixPropertyPlasmid() {
  if (plm_max) {
    memory->destroy(plm_x);
    memory->destroy(pre_x);
  }
}


/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPropertyPlasmid::grow_arrays(int nmax)
{
  FixProperty::grow_arrays(nmax);
  memory->grow(plm_x,nmax,plm_max*3,"fix_nufeb/property/plasmid:xpm");
  memory->grow(pre_x,nmax,3,"fix_nufeb/property/plasmid:xold");
}

/* ---------------------------------------------------------------------- */

int FixPropertyPlasmid::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ----------------------------------------------------------------------
   update cell age
------------------------------------------------------------------------- */

void FixPropertyPlasmid::biology_nufeb()
{
  compute();
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */

void FixPropertyPlasmid::compute()
{
  // record x
  for (int i = 0; i < atom->nlocal; i++) {
    pre_x[i][0] = atom->x[i][0];
    pre_x[i][1] = atom->x[i][1];
    pre_x[i][2] = atom->x[i][2];
  }
 // dump();
}

/* ----------------------------------------------------------------------
   update plasmid copy numbers in two daughter cells i, j
   called from fix_divide
------------------------------------------------------------------------- */
void FixPropertyPlasmid::update_arrays(int i, int j)
{
  if (atom->nlocal == atom->nmax) grow_arrays(atom->nmax);

  int ibac = atom->bacillus[i];
  double current_t = update->ntimestep * update->dt;
  int *dlist_plm, *dlist_fila, *tlist;

  memory->create(dlist_plm,plm_max,"fix nufeb/property/plasmid:dlist_plm");
  memory->create(dlist_fila,fila_max,"fix nufeb/property/plasmid:dlist_fila");
  memory->create(tlist,plm_max,"fix nufeb/property/plasmid:tlist");

  AtomVecBacillus::Bonus *ibouns = &avec->bonus[ibac];

  vprop[j] = 0;

  for (int m = 0; m < (int)vprop[i]; m++) {
    double px[3],pole1[3],pole2[3];
    double d,r;
    double idist = sqrt(2*atom->radius[i]*atom->radius[i])-(plm_dia*0.5);

    get_plasmid_coords(i, m, px, pre_x[i]);
    avec->get_pole_coords(i,pole1,pole2);
    distance_bt_pt_line(px,pole1,pole2,d);

    // plasmid transmission
    if (d > idist) {
      // move plasmid k to cell j
      int n = static_cast<int>(vprop[j]);
      dlist_plm[m] = 1;
      tlist[m] = n;
      if (rep_flag) nproteins[j][n] = nproteins[i][m];
      set_plm_x(j, n, px, atom->x[j]);
      vprop[j]++;
    } else {
      dlist_plm[m] = 0;
      set_plm_x(i, m, px, atom->x[i]);
    }
  }
  // delete broken filament
  if (par_flag) {
    nfilas[j] = 0;
    for (int f = 0; f < nfilas[i]; f++) {
      int njfila = nfilas[j];
      int m = fila[i][f][0];
      int n = fila[i][f][1];

      if (dlist_plm[m] != dlist_plm[n]) {
	dlist_fila[f] = 1;
      } else if (dlist_plm[m] && dlist_plm[n]) {
	dlist_fila[f] = 1;
	fila[j][njfila][0] = tlist[m];
	fila[j][njfila][1] = tlist[n];
	tfila[j][njfila] = tfila[i][f];
	nfilas[j]++;
      } else {
	dlist_fila[f] = 0;
      }
    }
    delete_filament(dlist_fila,i);
  }

  int k = 0;
  while (k < static_cast<int>(vprop[i])) {
    // remove plasmid k from i
    if (dlist_plm[k]) {
      copy_plasmid(i,k,static_cast<int>(vprop[i])-1,1);
      dlist_plm[k] = dlist_plm[static_cast<int>(vprop[i])-1];
      vprop[i]--;
    } else k++;
  }

  for (int m = 0; m < static_cast<int>(vprop[i]); m++) {
    relocate_plm_x(i,m);
  }

  for (int m = 0; m < static_cast<int>(vprop[j]); m++) {
    relocate_plm_x(j,m);
  }

  memory->destroy(dlist_plm);
  memory->destroy(dlist_fila);
  memory->destroy(tlist);
}

/* ----------------------------------------------------------------------
   delete filament
------------------------------------------------------------------------- */
void FixPropertyPlasmid::delete_filament(int *dflist, int i) {
  int k = 0;
  while (k < nfilas[i]) {
    int nfila = nfilas[i]-1;
    // remove plasmid k from i
    if (dflist[k]) {
      fila[i][k][0] = fila[i][nfila][0];
      fila[i][k][1] = fila[i][nfila][1];
      tfila[i][k] = tfila[i][nfila];

      dflist[k] = dflist[nfila];
      nfilas[i]--;
    } else k++;
  }
}

/* ----------------------------------------------------------------------
   get coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_plasmid_coords(int i, int j, double *px)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];

  double jplm_x[3];
  double v[3];
  double quat_pm[4];
  double v_length = bouns->length*0.5;

  jplm_x[0] = plm_x[i][j*3];
  jplm_x[1] = plm_x[i][j*3+1];
  jplm_x[2] = plm_x[i][j*3+2];

  v[0] = v[1] = v[2] = 0.0;

  if (jplm_x[0] > 0) {
    v[0] = v_length;
    get_quat(bouns->pole1,v,quat_pm);
  } else {
    v[0] = -v_length;
    get_quat(bouns->pole2,v,quat_pm);
  }

  double cpx[3];
  double p[3][3];

  MathExtra::quat_to_mat(quat_pm,p);
  MathExtra::matvec(p,jplm_x,cpx);
  MathExtra::quat_to_mat(bouns->quat,p);
  MathExtra::matvec(p,cpx,px);

  double *x = atom->x[bouns->ilocal];

  px[0] += x[0];
  px[1] += x[1];
  px[2] += x[2];
}

/* ----------------------------------------------------------------------
   get coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_plasmid_coords(int i, int j, double *xp, double *x)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];

  double jplm_x[3];
  double v[3];
  double quat_pm[4];
  double v_length = bouns->length*0.5;

  jplm_x[0] = plm_x[i][j*3];
  jplm_x[1] = plm_x[i][j*3+1];
  jplm_x[2] = plm_x[i][j*3+2];

  v[0] = v[1] = v[2] = 0.0;

  if (jplm_x[0] > 0) {
    v[0] = v_length;
    get_quat(bouns->pole1,v,quat_pm);
  } else {
    v[0] = -v_length;
    get_quat(bouns->pole2,v,quat_pm);
  }

  double cpx[3];
  double p[3][3];

  MathExtra::quat_to_mat(quat_pm,p);
  MathExtra::matvec(p,jplm_x,cpx);
  MathExtra::quat_to_mat(bouns->quat,p);
  MathExtra::matvec(p,cpx,xp);

  xp[0] += x[0];
  xp[1] += x[1];
  xp[2] += x[2];
}

/* ----------------------------------------------------------------------
   get quaternion from two vectors
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_quat(double *vec1, double *vec2, double *q){
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
   set relative coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::set_plm_x(int i, int j, double *xp, double *x)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];

  double jplm_x[3];
  double quat[4];
  double p[3][3];
  double cpx[3];

  quat[0] = bouns->quat[0];
  quat[1] = -bouns->quat[1];
  quat[2] = -bouns->quat[2];
  quat[3] = -bouns->quat[3];

  xp[0] -= x[0];
  xp[1] -= x[1];
  xp[2] -= x[2];

  MathExtra::quat_to_mat(quat,p);
  MathExtra::matvec(p,xp,cpx);

  double v[3];
  double quat_pm[4];
  double v_length = bouns->length*0.5;
  v[0] = v[1] = v[2] = 0.0;

  if (xp[0] > 0) {
    v[0] = v_length;
    get_quat(v,bouns->pole1,quat_pm);
  } else {
    v[0] = -v_length;
    get_quat(v,bouns->pole2,quat_pm);
  }

  MathExtra::quat_to_mat(quat_pm,p);
  MathExtra::matvec(p,cpx,jplm_x);

  plm_x[i][j*3] = jplm_x[0];
  plm_x[i][j*3+1] = jplm_x[1];
  plm_x[i][j*3+2] = jplm_x[2];
}

/* ----------------------------------------------------------------------
   get maximum/minimum coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_cell_boundary(double *xlimit, int i, int j)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
  double d = bouns->length/2;
  double r = bouns->diameter/2;

  // xhi
  xlimit[1] = d + r - plm_dia/2;
  // xlo
  xlimit[0] = -xlimit[1];
  // yhi zhi
  if (fabs(plm_x[i][j*3]) < d) xlimit[3] = xlimit[5] = r - plm_dia/2;
  else xlimit[3] = xlimit[5] = sqrt(r*r - (fabs(plm_x[i][j*3])-d)*(fabs(plm_x[i][j*3])-d)) - plm_dia/2;
  // yo zli
  xlimit[2] = -xlimit[3];
  xlimit[4] = -xlimit[5];
}

/* ----------------------------------------------------------------------
   reset position of plasmid j in bacillus i based on xlimit
------------------------------------------------------------------------- */
void FixPropertyPlasmid::relocate_plm_x(int i, int j) {
  double limit[6];
  int j0 = j*3;
  int j1 = j*3+1;
  int j2 = j*3+2;

  get_cell_boundary(limit, i, j);

  if (plm_x[i][j0] < limit[0])        // xlo
    plm_x[i][j0] = limit[0];
  else if (plm_x[i][j0] > limit[1])   // xhi
    plm_x[i][j0] = limit[1];
  if (plm_x[i][j1] < limit[2])        // ylo
    plm_x[i][j1] = limit[2];
  else if (plm_x[i][j1] > limit[3])   // yhi
    plm_x[i][j1] = limit[3];
  if (plm_x[i][j2] < limit[4])        // zlo
    plm_x[i][j2] = limit[4];
  else if (plm_x[i][j2] > limit[5])   // zhi
    plm_x[i][j2] = limit[5];
}

/* ----------------------------------------------------------------------
   distance between plasmid and host rod (line segment)
------------------------------------------------------------------------- */
void FixPropertyPlasmid::distance_bt_pt_line(double *q, double *xi1, double *xi2, double &d)
{
  double h[3], t;
  double vx = xi2[0] - xi1[0];
  double vy = xi2[1] - xi1[1];
  double vz = xi2[2] - xi1[2];

  double wx = q[0] - xi1[0];
  double wy = q[1] - xi1[1];
  double wz = q[2] - xi1[2];

  double c1 = dot(wx, wy, wz, vx, vy, vz);
  double c2 = dot(vx, vy, vz, vx, vy, vz);

  if (c1 <= 0) {
    t = 0;
    h[0] = xi1[0];
    h[1] = xi1[1];
    h[2] = xi1[2];

  } else if (c2 <= c1) {
    t = 1;
    h[0] = xi2[0];
    h[1] = xi2[1];
    h[2] = xi2[2];

  } else {
    t = c1 / c2;
    h[0] = xi1[0] + t * vx;
    h[1] = xi1[1] + t * vy;
    h[2] = xi1[2] + t * vz;
  }

  d = dist(q[0], q[1], q[2], h[0], h[1], h[2]);
}

/* ----------------------------------------------------------------------
   copy values within plasmid arrays
------------------------------------------------------------------------- */

void FixPropertyPlasmid::copy_plasmid(int i, int to, int from, int /*delflag*/)
{
  plm_x[i][to*3] = plm_x[i][from*3];
  plm_x[i][to*3+1] = plm_x[i][from*3+1];
  plm_x[i][to*3+2] = plm_x[i][from*3+2];

  if (rep_flag)
    nproteins[i][to] = nproteins[i][from];

  if (par_flag)
    for (int f = 0; f < nfilas[i]; f++) {
      if (fila[i][f][0] == from)
	fila[i][f][0] = to;
      if (fila[i][f][1] == from)
	fila[i][f][1] = to;
    }
}

/* ----------------------------------------------------------------------
   copy values within plasmid arrays
------------------------------------------------------------------------- */

void FixPropertyPlasmid::copy_arrays(int i, int j, int /*delflag*/)
{
  FixProperty::copy_arrays(i,j,0);

  for (int m = 0; m < plm_max*3; m++) {
    plm_x[j][m] = plm_x[i][m];
  }

  if (rep_flag) {
    for (int m = 0; m < plm_max; m++) {
      nproteins[j][m] = nproteins[i][m];
    }
  }

  if (par_flag) {
    for (int n = 0; n < fila_max; n++) {
      fila[j][n][0] = fila[i][n][0];
      fila[j][n][1] = fila[i][n][1];
      tfila[j][n] = tfila[i][n];
    }

    nfilas[j] = nfilas[i];
  }

  for (int m = 0; m < 3; m++) {
    pre_x[j][m] = pre_x[i][m];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyPlasmid::pack_exchange(int i, double *buf)
{
  int m = FixProperty::pack_exchange(i,buf);

  buf[m++] = vprop[i];

  for (int n = 0; n < plm_max*3; n++) {
    buf[m++] = plm_x[i][n];
  }

  if (rep_flag) {
    for (int n = 0; n < plm_max; n++) {
      buf[m++] = nproteins[i][n];
    }
  }

  if (par_flag) {
    for (int n = 0; n < fila_max; n++) {
      buf[m++] = fila[i][n][0];
      buf[m++] = fila[i][n][1];
      buf[m++] = tfila[i][n];
    }

    buf[m++] = nfilas[i];
  }

  buf[m++] = pre_x[i][0];
  buf[m++] = pre_x[i][1];
  buf[m++] = pre_x[i][2];

  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPropertyPlasmid::unpack_exchange(int nlocal, double *buf)
{
  int m = FixProperty::unpack_exchange(nlocal,buf);

  vprop[nlocal] = buf[m++];

  for (int n = 0; n < plm_max*3; n++) {
    plm_x[nlocal][n] = buf[m++];
  }

  if (rep_flag) {
    for (int n = 0; n < plm_max; n++) {
      nproteins[nlocal][n] = buf[m++];
    }
  }

  if (par_flag) {
    for (int n = 0; n < fila_max; n++) {
      fila[nlocal][n][0] = buf[m++];
      fila[nlocal][n][1] = buf[m++];
      tfila[nlocal][n] = buf[m++];
    }

    nfilas[nlocal] = buf[m++];
  }
  pre_x[nlocal][0] = buf[m++];
  pre_x[nlocal][1] = buf[m++];
  pre_x[nlocal][2] = buf[m++];

  return m;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropertyPlasmid::memory_usage()
{
  double bytes;

  bytes += atom->nmax*sizeof(double);
  bytes += atom->nmax*plm_max*sizeof(double)*3;

  return bytes;
}


/*------------------------------------------------------------------------- */
/*
void FixPropertyPlasmid::dump()
{
  if (update->ntimestep % 600) {
    myfile.open ("length7.txt", std::ios_base::app);
    for (int i = 0; i < atom->nlocal; i=i+10) {
      AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
      if (bouns->length > 6.98e-6 && bouns->length < 7.0e-6 && (int)vprop[i]) {
	for (int j = 0; j < (int)vprop[i]; j++) {
	    myfile << plm_x[i][j*3] << ", ";
	    myfile << plm_x[i][j*3+1] << ", ";
	    myfile << plm_x[i][j*3+2] << "";
	    myfile << "\n";
	}
      }
    }
    myfile.close();

    myfile.open ("length10.5.txt", std::ios_base::app);
    for (int i = 0; i < atom->nlocal; i=i+10) {
      AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
      if (bouns->length > 10.48e-6 && bouns->length < 10.5e-6 && (int)vprop[i]) {
	for (int j = 0; j < (int)vprop[i]; j++) {
	    myfile << plm_x[i][j*3] << ", ";
	    myfile << plm_x[i][j*3+1] << ", ";
	    myfile << plm_x[i][j*3+2] << "";
	    myfile << "\n";
	}
      }
    }
    myfile.close();

    myfile.open ("length14.txt", std::ios_base::app);
    for (int i = 0; i < atom->nlocal; i=i+10) {
      AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
      if (bouns->length > 13.95e-6 && bouns->length < 14e-6 && (int)vprop[i]) {
	for (int j = 0; j < (int)vprop[i]; j++) {
	    myfile << plm_x[i][j*3] << ", ";
	    myfile << plm_x[i][j*3+1] << ", ";
	    myfile << plm_x[i][j*3+2] << "";
	    myfile << "\n";
	}
      }
    }
    myfile.close();
  }

//  for bsub-wt-mm-30/n4-std, n4-nopar, n4-nopar-nonucle
    if (update->ntimestep % 490) {
      myfile.open ("length4.4.txt", std::ios_base::app);
      for (int i = 0; i < atom->nlocal; i=i+10) {
	AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
	if (bouns->length > 4.39e-6 && bouns->length < 4.3999e-6 && (int)vprop[i]) {
	  for (int j = 0; j < (int)vprop[i]; j++) {
	    double limit[6];
	    get_cell_boundary(limit, i, j);
	    if(plm_x[i][j*3] < limit[0] || plm_x[i][j*3] > limit[1] ||
	       plm_x[i][j*3+1] < limit[2] || plm_x[i][j*3+1] > limit[3] ||
	       plm_x[i][j*3+2] < limit[4] || plm_x[i][j*3+2] > limit[5]){
		continue;
	    }
	    myfile << plm_x[i][j*3] << ", ";
	    myfile << plm_x[i][j*3+1] << ", ";
	    myfile << plm_x[i][j*3+2] << "";
	    myfile << "\n";
	  }
	}
      }
      myfile.close();
    }


//  for bsub-wt-mm-30/single-1plm and single-2plm
    myfile.open ("trajectory.txt", std::ios_base::app);
    for (int i = 0; i < atom->nlocal; i++) {
      // initialise plasmid position
      for (int j = 0; j < (int)vprop[i]; j++) {
  	myfile << plm_x[i][j*3] << ", ";
      }
      myfile << "\n";
    }
    myfile.close();
}*/
