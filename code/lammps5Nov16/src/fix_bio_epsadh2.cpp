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

#include "fix_bio_epsadh2.h"

#include <math.h>
#include <string.h>

#include "atom.h"
#include "atom_vec_bio.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pointers.h"
#include "variable.h"


using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 10000;
/* ---------------------------------------------------------------------- */

FixEPSAdh2::FixEPSAdh2(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 6) error->all(FLERR,"Illegal fix eps adhesion command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix eps adhesion command: calling steps should be positive integer");
  flag = force->inumeric(FLERR, arg[5]);
  if (flag != 1 && flag != 2) error->all(FLERR,"Illegal fix eps adhesion command: undefined model");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var,&arg[4][2]);
}

FixEPSAdh2::~FixEPSAdh2()
{
  delete [] var;
}

/* ---------------------------------------------------------------------- */

int FixEPSAdh2::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEPSAdh2::init()
{
  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR,"Variable name for fix eps adhesion does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR,"Variable for fix eps adhesion is invalid style");
}

/* ---------------------------------------------------------------------- */

void FixEPSAdh2::post_force(int vflag)
{
  double ke = input->variable->compute_equal(ivar);

  int i,ii=0,j=0,jj=0,*numneigh,inum,**firstneigh;
  int jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double outerRadi,outerRadj,radsum,rsq,r, rinv;
  double ccel,ccelx ,ccely,ccelz;
  int *jlist,*ilist;
  double del;
  double **f = atom->f;
  double **x = atom->x;
  double *outerRadius = avec->outerRadius;
  double *outerMass = avec->outerMass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int *mask = atom->mask;
  double epsMassi, epsMassj, massSum;

  // loop over neighbors of my atoms
  
  for (i = 0; i < nlocal; i++){
	  if (!(mask[i] & groupbit)) continue;
	  xtmp = x[i][0];
	  ytmp = x[i][1];
	  ztmp = x[i][2];

	  outerRadi = outerRadius[i];
	  jlist = firstneigh[i];
	  jnum = numneigh[i];

    if (atom->mask[i] == avec->maskEPS) epsMassi = rmass[i];
    else epsMassi = outerMass[i];

    for(int j = 0; j < atom->nlocal; j++){
      if(i != j){
        double xd = atom->x[i][0] - atom->x[j][0];
        double yd = atom->x[i][1] - atom->x[j][1];
        double zd = atom->x[i][2] - atom->x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (atom->radius[i] + atom->radius[j] + 1e-6) * (atom->radius[i] + atom->radius[j]+ 1e-6);
        if (rsq <= cut) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;

          outerRadj = outerRadius[j];

          if (atom->mask[j] == avec->maskEPS) epsMassj = rmass[j];
          else epsMassj = outerMass[j];

          radsum = outerRadi + outerRadj;
          massSum = epsMassi + epsMassj;

          if(flag == 1){
            if (rsq < 4 * radsum * radsum) {
              r = sqrt(rsq);
              del = r - radsum;
              rinv = 1/r;
              ccel = -massSum*ke*del;
            }
          }else if(flag == 2){
            if ((rsq < 4 * radsum * ke) && (rsq > radsum)){
              r = sqrt(rsq);
              del = r - radsum;
              rinv = 1/r;
              ccel = -massSum*ke*(radsum/r)*(radsum/r);
            }
          }

          ccelx = delx*ccel*rinv ;
          ccely = dely*ccel*rinv ;
          ccelz = delz*ccel*rinv ;

          f[i][0] += ccelx;
          f[i][1] += ccely;
          f[i][2] += ccelz;

          if (newton_pair || j < nlocal) {
            f[j][0] -= ccelx;
            f[j][1] -= ccely;
            f[j][2] -= ccelz;
          }
        }
      }
	  }
  }
}
