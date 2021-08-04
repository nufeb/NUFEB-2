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

#include "fix_property_plasmid.h"

#include <cstdlib>
#include <cstring>
#include <math.h>
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "grid.h"
#include "atom_vec_bacillus.h"
#include "update.h"
#include "math_const.h"
#include "math_extra.h"
#include "random_park.h"
#include <random>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define FILAMENT_VEL 0.026e-6
#define CUTOFF 3e-7
#define DELTA 1.005

/* ---------------------------------------------------------------------- */

FixPropertyPlasmid::FixPropertyPlasmid(LAMMPS *lmp, int narg, char **arg) :
  FixProperty(lmp, narg, arg),
  xpm(nullptr), fila(nullptr), pre_x(nullptr), nproteins(nullptr),
  tfila(nullptr), nfilas(nullptr)
{
  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"fix nufeb/property/plasmid requires "
      "atom style bacillus");

  if (narg < 3) error->all(FLERR,"Illegal fix nufeb/property/plasmid command");

  compute_flag = 1;
  scalar_flag = 1;

  mean_protein = 20;
  init_protein = 0;
  pmax = 5;
  fmax = 0;
  pinit = 0;
  dia = 5e-7;
  diff_coef = 4e-15;
  dt = 1;
  alpha = 1.0;
  ftime = 60;

  seed = utils::inumeric(FLERR,arg[3],true,lmp);
  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "replicate") == 0) {
      mean_protein = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      init_protein = utils::inumeric(FLERR,arg[iarg+2],true,lmp);
      alpha = utils::numeric(FLERR,arg[iarg+3],true,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "max") == 0) {
      pmax = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "init") == 0) {
      pinit = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dia") == 0) {
      dia = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "diffusion") == 0) {
      diff_coef = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dt") == 0) {
      dt = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "ftime") == 0) {
      ftime = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/property/plasmid command");
    }
  }

  if (pinit > pmax) error->all(FLERR,"Illegal fix nufeb/property/plasmid command: "
      "initial plasmid cannot be more than maximum plasmid number");

  size_peratom_cols = 0;
  fmax = pmax * (pmax - 1) / 2;
  grow_arrays(atom->nmax);

  std::random_device rd;
  std::mt19937 eng(rd());
  for (int i = 0; i < atom->nlocal; i++) {
    AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
    vprop[i] = pinit;
    // initialise plasmid position
    for (int j = 0; j < (int)vprop[i]; j++) {
      double limit[6];

      int jx0 = j*3;
      int jx1 = j*3+1;
      int jx2 = j*3+2;

      // assign x
      limit[1] = (bouns->length + bouns->diameter - dia) / 2;
      limit[0] = -limit[1];

      std::uniform_real_distribution<double> x(limit[0], limit[1]);
      xpm[i][jx0] = x(eng);

      get_xlimit(limit, i, j);

      // assign y, z
      std::uniform_real_distribution<double> y(limit[2], limit[3]);
      std::uniform_real_distribution<double> z(limit[4], limit[5]);
      xpm[i][jx1] = y(eng);
      xpm[i][jx2] = z(eng);

    //fixed positions: 2 plms
//      if (j == 0) xpm[i][jx0] = 0.75e-6;
//      if (j == 1) xpm[i][jx0] = -0.75e-6;

      //fixed positions: 3 plms
//      if (j == 0) xpm[i][jx0] = 1.3e-6;
//      if (j == 1) xpm[i][jx0] = -1.2e-6;
//      if (j == 2) xpm[i][jx0] = 0.45e-6;

      nproteins[i][j] = init_protein;
    }

    for (int f = 0; f < fmax; f++) {
      tfila[i][f] = 0.0;
      fila[i][f][0] = -1;
      fila[i][f][1] = -1;
    }

    nfilas[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */
FixPropertyPlasmid:: ~FixPropertyPlasmid() {
  memory->destroy(nfilas);

  if (pmax) {
    memory->destroy(xpm);
    memory->destroy(pre_x);
    memory->destroy(nproteins);
  }

  if (fmax) {
    memory->destroy(fila);
    memory->destroy(tfila);
  }

  delete random;
}


/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPropertyPlasmid::grow_arrays(int nmax)
{
  FixProperty::grow_arrays(nmax);
  memory->grow(xpm,nmax,pmax*3,"fix_nufeb/property/plasmid:xpm");
  memory->grow(pre_x,nmax,3,"fix_nufeb/property/plasmid:xold");
  memory->grow(nproteins,nmax,pmax,"fix_nufeb/property/plasmid:ninitiator");

  memory->grow(fila,nmax,fmax,2,"fix_nufeb/property/plasmid:fila");
  memory->grow(nfilas,nmax,"fix_nufeb/property/plasmid:nfilas");
  memory->grow(tfila,nmax,fmax,"fix_nufeb/property/plasmid:tfila");
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */

void FixPropertyPlasmid::compute()
{
  double t = 0;
  while (t < update->dt) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit)
	motility(i);
    }
    t+=dt;
  }

  for (int i = 0; i < atom->nlocal; i++){
    if (atom->mask[i] & groupbit) {
      replication(i);
    }
  }

  // record x
  for (int i = 0; i < atom->nlocal; i++) {
    pre_x[i][0] = atom->x[i][0];
    pre_x[i][1] = atom->x[i][1];
    pre_x[i][2] = atom->x[i][2];
  }

//  myfile.open ("trace.txt", std::ios_base::app);
//  for (int i = 0; i < atom->nlocal; i++) {
//    // initialise plasmid position
//    for (int j = 0; j < (int)vprop[i]; j++) {
//	myfile << xpm[i][j*3] << ", ";
//    }
//    myfile << "\n";
//  }
//  myfile.close();
}

/* ----------------------------------------------------------------------
   update plasmid position in atom i
------------------------------------------------------------------------- */

void FixPropertyPlasmid::motility(int i) {
  int *dflist;

  memory->create(dflist,fmax,"fix nufeb/property/plasmid:dlist");


  for (int m = 0; m < (int)vprop[i]; m++) {
    double pos[3];
    double xlimit[3];

    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;

    // create filament between plasmids n and m if collision
    for (int n = 0; n < (int)vprop[i]; n++) {
      if (!ftime) break;

      int nfila = nfilas[i];
      if (m == n) continue;
      int skip = 0;
      for (int f = 0; f < nfila; f++) {
        if ((fila[i][f][0] == m && fila[i][f][1] == n) ||
            (fila[i][f][0] == n && fila[i][f][1] == m)) {
          skip = 1;
          break;
        }
      }
      if (skip) continue;

      int n0 = n*3;
      int n1 = n*3+1;
      int n2 = n*3+2;

      double rsq;
      rsq = (xpm[i][n2]-xpm[i][m2])*(xpm[i][n2]-xpm[i][m2]) +
	  (xpm[i][n1]-xpm[i][m1])*(xpm[i][n1]-xpm[i][m1]) +
	  (xpm[i][n0]-xpm[i][m0])*(xpm[i][n0]-xpm[i][m0]);

      if (rsq < dia*dia + CUTOFF*CUTOFF) {
	if (xpm[i][m0] > xpm[i][n0]) {
	  fila[i][nfila][0] = n;
	  fila[i][nfila][1] = m;
	} else {
	  fila[i][nfila][0] = m;
	  fila[i][nfila][1] = n;
	}

	tfila[i][nfila] = ftime;
	nfilas[i]++;
      }
    }

    int dir = 0;
    int con = 0;

    // update filament attribute
    for (int f = 0; f < nfilas[i]; f++) {
      int f0 = fila[i][f][0];
      int f1 = fila[i][f][1];

      if (f0 == m || f1 == m) {
	con = 1;
	if (tfila[i][f] > 0) {
	  dflist[f] = 0;
	  tfila[i][f] -= dt/2;

	  if (f0 == m) dir--;
	  else dir++;
	} else {
	  dflist[f] = 1;
	}
      }
    }

    // pushing by ParM filament
    if (con) {
      if (dir < 0) {
	pos[0] = xpm[i][m0] - FILAMENT_VEL * dt;
	pos[1] = xpm[i][m1];
	pos[2] = xpm[i][m2];
      } else if (dir > 0){
	pos[0] = xpm[i][m0] + FILAMENT_VEL * dt;
	pos[1] = xpm[i][m1];
	pos[2] = xpm[i][m2];
      } else {
	pos[0] = xpm[i][m0];
	pos[1] = xpm[i][m1];
	pos[2] = xpm[i][m2];
      }
    } else {
    // Brownian motion
      pos[0] = xpm[i][m0] + sqrt(2*diff_coef*dt)*random->gaussian();
      pos[1] = xpm[i][m1] + sqrt(2*diff_coef*dt)*random->gaussian();
      pos[2] = xpm[i][m2] + sqrt(2*diff_coef*dt)*random->gaussian();
    }

    xpm[i][m0] = pos[0];
    xpm[i][m1] = pos[1];
    xpm[i][m2] = pos[2];

    relocate_limit(i,m);
  }

  delete_filament(dflist,i);
  memory->destroy(dflist);
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
   update plasmid copy number
------------------------------------------------------------------------- */
void FixPropertyPlasmid::replication(int i) {
  double current_t = update->ntimestep * update->dt;
  double next_t = (update->ntimestep + 1) * update->dt;
  double *ixpm = xpm[i];
  const int cell = grid->cell(atom->x[i]);
  double growth = grid->growth[igroup][cell][0];

  // update initiator proteins
  for (int m = 0; m < (int)vprop[i]; m++)
    nproteins[i][m] += alpha * growth * atom->rmass[i] * update->dt;

  for (int m = 0; m < (int)vprop[i]; m++) {
    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;
    double max = mean_protein + mean_protein*0.1*random->gaussian();

    if (nproteins[i][m] > max && vprop[i] < pmax) {
      double ilimit[3];

      int n = vprop[i];
      int n0 = n*3;
      int n1 = n*3+1;
      int n2 = n*3+2;

      double theta = random->uniform() * 2 * MY_PI;
      double phi = random->uniform() * (MY_PI);

      ixpm[n0] = ixpm[m0] + (dia * cos(theta) * sin(phi) * DELTA);
      ixpm[n1] = ixpm[m1] + (dia * sin(theta) * sin(phi) * DELTA);
      ixpm[n2] = ixpm[m2] + (dia * cos(phi) * DELTA);

      relocate_limit(i,m);
      vprop[i]++;

      nproteins[i][m] = 0.0;
      nproteins[i][n] = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   get coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_plasmid_coords(int i, int j, double *xp)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];

  double jxpm[3];
  double v[3];
  double quat_pm[4];
  double v_length = bouns->length*0.5;

  jxpm[0] = xpm[i][j*3];
  jxpm[1] = xpm[i][j*3+1];
  jxpm[2] = xpm[i][j*3+2];

  v[0] = v[1] = v[2] = 0.0;

  if (jxpm[0] > 0) {
    v[0] = v_length;
    get_quat(bouns->pole1,v,quat_pm);
  } else {
    v[0] = -v_length;
    get_quat(bouns->pole2,v,quat_pm);
  }

  double cpx[3];
  double p[3][3];

  MathExtra::quat_to_mat(quat_pm,p);
  MathExtra::matvec(p,jxpm,cpx);
  MathExtra::quat_to_mat(bouns->quat,p);
  MathExtra::matvec(p,cpx,xp);

  double *x = atom->x[bouns->ilocal];

  xp[0] += x[0];
  xp[1] += x[1];
  xp[2] += x[2];
}

/* ----------------------------------------------------------------------
   get coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_plasmid_coords(int i, int j, double *xp, double *x)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];

  double jxpm[3];
  double v[3];
  double quat_pm[4];
  double v_length = bouns->length*0.5;

  jxpm[0] = xpm[i][j*3];
  jxpm[1] = xpm[i][j*3+1];
  jxpm[2] = xpm[i][j*3+2];

  v[0] = v[1] = v[2] = 0.0;

  if (jxpm[0] > 0) {
    v[0] = v_length;
    get_quat(bouns->pole1,v,quat_pm);
  } else {
    v[0] = -v_length;
    get_quat(bouns->pole2,v,quat_pm);
  }

  double cpx[3];
  double p[3][3];

  MathExtra::quat_to_mat(quat_pm,p);
  MathExtra::matvec(p,jxpm,cpx);
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
void FixPropertyPlasmid::set_plasmid_xpm(int i, int j, double *xp, double *x)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];

  double jxpm[3];
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
  MathExtra::matvec(p,cpx,jxpm);

  xpm[i][j*3] = jxpm[0];
  xpm[i][j*3+1] = jxpm[1];
  xpm[i][j*3+2] = jxpm[2];
}

/* ----------------------------------------------------------------------
   get maximum/minimum coordinate for plasmid j in bacillus i
------------------------------------------------------------------------- */
void FixPropertyPlasmid::get_xlimit(double *xlimit, int i, int j)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
  double d = bouns->length/2;
  double r = bouns->diameter/2;

  // xhi
  xlimit[1] = d + r - dia/2;
  // xlo
  xlimit[0] = -xlimit[1];
  // yhi zhi
  if (fabs(xpm[i][j*3]) < d) xlimit[3] = xlimit[5] = r - dia/2;
  else xlimit[3] = xlimit[5] = sqrt(r*r - (fabs(xpm[i][j*3])-d)*(fabs(xpm[i][j*3])-d)) - dia/2;
  // yo zli
  xlimit[2] = -xlimit[3];
  xlimit[4] = -xlimit[5];
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
  int *dlist, *dflist, *tlist;

  memory->create(dlist,pmax,"fix nufeb/property/plasmid:dlist");
  memory->create(dflist,fmax,"fix nufeb/property/plasmid:dflist");
  memory->create(tlist,pmax,"fix nufeb/property/plasmid:jlist");

  AtomVecBacillus::Bonus *ibouns = &avec->bonus[ibac];

  vprop[j] = 0;
  nfilas[j] = 0;

  for (int m = 0; m < (int)vprop[i]; m++) {
    double xp[3],xp1[3],xp2[3];
    double d,r;
    double idist = sqrt(2*atom->radius[i]*atom->radius[i])-(dia*0.5);

    get_plasmid_coords(i, m, xp, pre_x[i]);
    avec->get_pole_coords(i,xp1,xp2);
    distance_bt_pt_line(xp,xp1,xp2,d);

    // plasmid transmission
    if (d > idist) {
      // move plasmid k to cell j
      int n = (int)vprop[j];
      dlist[m] = 1;
      tlist[m] = n;
      nproteins[j][n] = nproteins[i][m];
      set_plasmid_xpm(j, n, xp, atom->x[j]);
      vprop[j]++;
    } else {
      dlist[m] = 0;
      set_plasmid_xpm(i, m, xp, atom->x[i]);
    }
  }

  // delete broken filament
  for (int f = 0; f < nfilas[i]; f++) {
    int njfila = nfilas[j];
    int m = fila[i][f][0];
    int n = fila[i][f][1];

    if (dlist[m] != dlist[n]) {
      dflist[f] = 1;
    } else if (dlist[m] && dlist[n]) {
      dflist[f] = 1;
      fila[j][njfila][0] = tlist[m];
      fila[j][njfila][1] = tlist[n];
      tfila[j][njfila] = tfila[i][f];
      nfilas[j]++;
    } else {
      dflist[f] = 0;
    }
  }

  delete_filament(dflist,i);

  int k = 0;
  while (k < (int)vprop[i]) {
    // remove plasmid k from i
    if (dlist[k]) {
      copy_plasmid(i,k,(int)vprop[i]-1,1);
      dlist[k] = dlist[(int)vprop[i]-1];
      vprop[i]--;
    } else k++;
  }

  for (int m = 0; m < (int)vprop[i]; m++) {
    relocate_limit(i,m);
  }

  for (int m = 0; m < (int)vprop[j]; m++) {
    relocate_limit(j,m);
  }

  memory->destroy(dlist);
  memory->destroy(dflist);
  memory->destroy(tlist);
}

/* ----------------------------------------------------------------------
   reset position of plasmid j in bacillus i based on xlimit
------------------------------------------------------------------------- */
void FixPropertyPlasmid::relocate_limit(int i, int j) {
  double limit[6];
  int j0 = j*3;
  int j1 = j*3+1;
  int j2 = j*3+2;

  get_xlimit(limit, i, j);

  if (xpm[i][j0] < limit[0])        // xlo
    xpm[i][j0] = limit[0];
  else if (xpm[i][j0] > limit[1])   // xhi
    xpm[i][j0] = limit[1];
  if (xpm[i][j1] < limit[2])        // ylo
    xpm[i][j1] = limit[2];
  else if (xpm[i][j1] > limit[3])   // yhi
    xpm[i][j1] = limit[3];
  if (xpm[i][j2] < limit[4])        // zlo
    xpm[i][j2] = limit[4];
  else if (xpm[i][j2] > limit[5])   // zhi
    xpm[i][j2] = limit[5];
}

/* ----------------------------------------------------------------------
   distance between plasmid and
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
   compute average plasmid copy number
------------------------------------------------------------------------- */
double FixPropertyPlasmid::compute_scalar() {
  double result = 0.0;
  int n = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((atom->mask[i] & groupbit)) {
      result += vprop[i];
      n++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_INT, MPI_SUM, world);

  if (n > 0) result /= n;
  else result = 0;

  return result;
}

/* ----------------------------------------------------------------------
   copy values within plasmid arrays
------------------------------------------------------------------------- */

void FixPropertyPlasmid::copy_plasmid(int i, int to, int from, int /*delflag*/)
{
  xpm[i][to*3] = xpm[i][from*3];
  xpm[i][to*3+1] = xpm[i][from*3+1];
  xpm[i][to*3+2] = xpm[i][from*3+2];

  nproteins[i][to] = nproteins[i][from];

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

  for (int m = 0; m < pmax*3; m++) {
    xpm[j][m] = xpm[i][m];
  }

  for (int m = 0; m < pmax; m++) {
    nproteins[j][m] = nproteins[i][m];
  }

  for (int n = 0; n < fmax; n++) {
    fila[j][n][0] = fila[i][n][0];
    fila[j][n][1] = fila[i][n][1];
    tfila[j][n] = tfila[i][n];
  }

  nfilas[j] = nfilas[i];

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

  for (int n = 0; n < pmax*3; n++) {
    buf[m++] = xpm[i][n];
  }

  for (int n = 0; n < pmax; n++) {
    buf[m++] = nproteins[i][n];
  }

  for (int n = 0; n < fmax; n++) {
    buf[m++] = fila[i][n][0];
    buf[m++] = fila[i][n][1];
    buf[m++] = tfila[i][n];
  }

  buf[m++] = nfilas[i];
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

  for (int n = 0; n < pmax*3; n++) {
    xpm[nlocal][n] = buf[m++];
  }

  for (int n = 0; n < pmax; n++) {
    nproteins[nlocal][n] = buf[m++];
  }

  for (int n = 0; n < fmax; n++) {
    fila[nlocal][n][0] = buf[m++];
    fila[nlocal][n][1] = buf[m++];
    tfila[nlocal][n] = buf[m++];
  }

  nfilas[nlocal] = buf[m++];
  pre_x[nlocal][0] = buf[m++];
  pre_x[nlocal][1] = buf[m++];
  pre_x[nlocal][2] = buf[m++];

  return m;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixPropertyPlasmid::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = pmax*3 + 1;
  buf[m++] = vprop[i];

  for (int n = 0; n < pmax*3; n++)
    buf[m++] = xpm[i][n];

  return pmax*3 + 1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixPropertyPlasmid::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  vprop[nlocal] = extra[nlocal][m++];
  for (int n = 0; n < pmax*3; n++)
    xpm[nlocal][n] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixPropertyPlasmid::maxsize_restart()
{
  return pmax*3 + 1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixPropertyPlasmid::size_restart(int /*nlocal*/)
{
  return pmax*3 + 1;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPropertyPlasmid::memory_usage()
{
  double bytes;
  bytes += atom->nmax*sizeof(double);
  bytes += atom->nmax*pmax*sizeof(double);

  return bytes;
}
