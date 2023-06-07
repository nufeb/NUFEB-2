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
#include "grid.h"
#include "atom_masks.h"

#include "fix_plasmid_partition.h"
#include "fix_property_plasmid.h"
#include "fix_divide.h"
#include "fix_divide_bacillus.h"
#include "fix_divide_bacillus_minicell.h"

#include "fix.h"
#include "random_park.h"
#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define CUTOFF 0.5e-6
#define NUCLEOID_DIA_RATIO 0.65
#define NUCLEOID_LEN_RATIO 0.65
#define BETA 1.01

/* ---------------------------------------------------------------------- */

FixPlasmidPartition::FixPlasmidPartition(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR, "Illegal nufeb/plasmid/replication command");

  fix_plm = nullptr;
  nfilas = nullptr;
  tfila = nullptr;
  fila = nullptr;

  diff_coef = 4e-15;
  dt = 1;
  tmax_fila = 60;
  v_fila = 0.026e-6;
  nucleoid_flag = 1;
  fila_max = 0;
  divflag = -1;

  seed = utils::inumeric(FLERR,arg[3],true,lmp);
  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);

  avec = (AtomVecBacillus *) atom->style_match("bacillus");
  if (!avec) error->all(FLERR,"fix nufeb/property/plasmid requires "
      "atom style bacillus");

  auto fixlist = modify->get_fix_by_style("^nufeb/property/plasmid");
  if (fixlist.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/property/plasmid defined for fix plasmid/replicatio");
  fix_plm = dynamic_cast<FixPropertyPlasmid *>(fixlist.front());

  auto fixlist2 = modify->get_fix_by_style("^nufeb/division/bacillus/minicell");
  auto fixlist3 = modify->get_fix_by_style("^nufeb/division/bacillus");

  if (fixlist2.size() != 1 || fixlist3.size() != 1)
    error->all(FLERR, "There must be exactly one fix nufeb/division/bacillus/* defined for fix plasmid/partition");
  if (fixlist2.size() == 1) {
    divflag = 0;
    fix_div_mini = static_cast<FixDivideBacillusMinicell *>(fixlist2.front());
  } else if (fixlist3.size() == 1) {
    fix_div = static_cast<FixDivideBacillus *>(fixlist3.front());
    divflag = 1;
  }

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "diffusion") == 0) {
      diff_coef = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "dt") == 0) {
      dt = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "fdur") == 0) {
      tmax_fila = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }  else if (strcmp(arg[iarg], "fvel") == 0) {
      v_fila = utils::numeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    }  else if (strcmp(arg[iarg], "nflag") == 0) {
      nucleoid_flag = utils::inumeric(FLERR,arg[iarg+1],true,lmp);
      iarg += 2;
    } else {
      error->all(FLERR,"Illegal fix nufeb/plasmid/replication command");
    }
  }

  fix_plm->par_flag = 1;
}

/* ---------------------------------------------------------------------- */
FixPlasmidPartition:: ~FixPlasmidPartition() {
  memory->destroy(nfilas);
  memory->destroy(fila);
  memory->destroy(tfila);

  delete random;
}

/* ---------------------------------------------------------------------- */
void FixPlasmidPartition::init() {
  fila_max = fix_plm->fila_max;
  grow_arrays(atom->nmax);

  for (int i = 0; i < atom->nlocal; i++) {
    // initialise protein number
    for (int j = 0; j < static_cast<int>(fix_plm->vprop[i]); j++) {
      for (int f = 0; f < fix_plm->fila_max; f++) {
	tfila[i][f] = 0.0;
	fila[i][f][0] = -1;
	fila[i][f][1] = -1;
      }
    }
    nfilas[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixPlasmidPartition::grow_arrays(int nmax)
{
  memory->grow(fila,nmax,fila_max,2,"fix_nufeb/property/plasmid:fila");
  memory->grow(nfilas,nmax,"fix_nufeb/property/plasmid:nfilas");
  memory->grow(tfila,nmax,fila_max,"fix_nufeb/property/plasmid:tfila");

  fix_plm->nfilas = nfilas;
  fix_plm->fila = fila;
  fix_plm->tfila = tfila;
}


/* ---------------------------------------------------------------------- */

int FixPlasmidPartition::setmask()
{
  int mask = 0;
  mask |= BIOLOGY_NUFEB;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPlasmidPartition::biology_nufeb()
{
  compute();
}

/* ---------------------------------------------------------------------- */

void FixPlasmidPartition::compute()
{
  double t = 0;
  while (t < update->dt) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (atom->mask[i] & groupbit)
	partition(i);
    }
    t+=dt;
  }
}

/* ----------------------------------------------------------------------
   update plasmid copy number
------------------------------------------------------------------------- */
void FixPlasmidPartition::partition(int i)
{
  int *dflist;

  memory->create(dflist,fila_max,"fix nufeb/property/plasmid:dlist");
  double **plm_x = fix_plm->plm_x;
  double plm_dia = fix_plm->plm_dia;

  for (int m = 0; m < static_cast<int>(fix_plm->vprop[i]); m++) {
    double pos[3];
    double xlimit[3];

    int m0 = m*3;
    int m1 = m*3+1;
    int m2 = m*3+2;

    // create filament between plasmids n and m if collision
    for (int n = 0; n < static_cast<int>(fix_plm->vprop[i]); n++) {
      if (!tmax_fila) break;

      int nfila = nfilas[i];
      if (m == n) continue;
      int skip = 0;
      for (int f = 0; f < nfila; f++) {
        if ((fila[i][f][0] == m && fila[i][f][1] == n) ||
            (fila[i][f][0] == n && fila[i][f][1] == m)) {
          skip = 1;
          break;
        }
        // one filament per plasmid
        if (fila[i][f][0] == m || fila[i][f][1] == m ||
            fila[i][f][0] == n || fila[i][f][1] == n) {
	  skip = 1;
          break;
        }
      }
      if (skip) continue;

      int n0 = n*3;
      int n1 = n*3+1;
      int n2 = n*3+2;

      double rsq;
      rsq = (plm_x[i][n2]-plm_x[i][m2])*(plm_x[i][n2]-plm_x[i][m2]) +
	  (plm_x[i][n1]-plm_x[i][m1])*(plm_x[i][n1]-plm_x[i][m1]) +
	  (plm_x[i][n0]-plm_x[i][m0])*(plm_x[i][n0]-plm_x[i][m0]);

      if (rsq < plm_dia*plm_dia + CUTOFF*CUTOFF) {
	if (plm_x[i][m0] > plm_x[i][n0]) {
	  fila[i][nfila][0] = n;
	  fila[i][nfila][1] = m;
	} else {
	  fila[i][nfila][0] = m;
	  fila[i][nfila][1] = n;
	}

	tfila[i][nfila] = tmax_fila;
	nfilas[i]++;
      }
    }

    // for multilinked plasmid
    int orient = 0;
    int link = 0;

    // update filament attribute
    for (int f = 0; f < nfilas[i]; f++) {
      int f0 = fila[i][f][0];
      int f1 = fila[i][f][1];

      if (f0 == m || f1 == m) {
	link = 1;
	if (tfila[i][f] > 0) {
	  dflist[f] = 0;
	  tfila[i][f] -= dt/2;

	  if (f0 == m) orient--;
	  else orient++;
	} else {
	  dflist[f] = 1;
	}
      }
    }

    // pushing by ParM filament
    if (link) {
      if (orient < 0) {
	pos[0] = plm_x[i][m0] - v_fila * dt;
	pos[1] = plm_x[i][m1];
	pos[2] = plm_x[i][m2];
      } else if (orient > 0){
	pos[0] = plm_x[i][m0] + v_fila * dt;
	pos[1] = plm_x[i][m1];
	pos[2] = plm_x[i][m2];
      } else {
	pos[0] = plm_x[i][m0];
	pos[1] = plm_x[i][m1];
	pos[2] = plm_x[i][m2];
      }
    } else {
      double dcx, dcy, dcz;
      int nucl;
      dcx = dcy = dcz = diff_coef;
      // check if plasmid xpm0 is in nucleoid area
      if (nucleoid_flag ) {
	nucl = check_nucleoid(i, m, plm_x[i][m0]);
	if (nucl) {
	  dcx = diff_coef;
	  dcy = diff_coef * 0.1;
	  dcz = diff_coef * 0.1;
	}
      }

      // Brownian motion
      pos[0] = plm_x[i][m0] + sqrt(2*dcx*dt)*random->gaussian();
      pos[1] = plm_x[i][m1] + sqrt(2*dcy*dt)*random->gaussian();
      pos[2] = plm_x[i][m2] + sqrt(2*dcz*dt)*random->gaussian();

      if (nucleoid_flag) {
	int nucl1 = check_nucleoid(i, m, pos[0]);
	double rsq = atom->radius[i] * atom->radius[i] * NUCLEOID_DIA_RATIO;

	// check if plasmid is entering nucleoid
	if (!nucl && nucl1 && (pos[1]*pos[1] + pos[2]*pos[2]) < rsq) {
	  pos[0] = plm_x[i][m0];
	}
      }
    }

    plm_x[i][m0] = pos[0];
    plm_x[i][m1] = pos[1];
    plm_x[i][m2] = pos[2];

    fix_plm->relocate_plm_x(i,m);
  }

  fix_plm->delete_filament(dflist,i);
  memory->destroy(dflist);
}

/* ----------------------------------------------------------------------
   check if plasmid j is in nucleoid area
------------------------------------------------------------------------- */
int FixPlasmidPartition::check_nucleoid(int i, int j, double xpm0)
{
  if (atom->bacillus[i] < 0)
    error->one(FLERR,"Assigning bacillus parameters to non-bacillus atom");

  AtomVecBacillus::Bonus *bouns = &avec->bonus[atom->bacillus[i]];
  double lb;
  if (divflag)
    lb = fix_div_mini->maxlength;
  else
    lb = fix_div->maxlength;
  double l = bouns->length;
  double ln = lb * NUCLEOID_LEN_RATIO;

  if (l < 2*lb*BETA) {
    double n1_lo = -0.5 * (l-lb+ln);
    double n1_hi = -0.5 * (l-lb-ln);

    if ((xpm0 > n1_lo && xpm0 < n1_hi) || (xpm0 > -n1_hi && xpm0 < -n1_lo))
      return 1;
    else
      return 0;
  } else {
    double n1_lo = -0.5 * (0.5*l+ln);
    double n1_hi = -0.5 * (0.5*l-ln);
    double n2_lo = -0.5 * (2*lb-0.5*l+ln);
    double n2_hi = -0.5 * (2*lb-0.5*l-ln);

    if ((xpm0 > n1_lo && xpm0 < n1_hi) || (xpm0 > n2_lo && xpm0 < n2_hi) ||
	(xpm0 > -n1_hi && xpm0 < -n1_lo) || (xpm0 > -n2_hi && xpm0 < -n2_lo))
      return 1;
    else
      return 0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixPlasmidPartition::memory_usage()
{
  double bytes;

  bytes += atom->nmax*sizeof(double);
  bytes += atom->nmax*fila_max*sizeof(double)*3;

  return bytes;
}
