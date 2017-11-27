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

#include "fix_bio_eps_extract.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "atom.h"
#include "atom_vec_bio.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lmptype.h"
#include "math_const.h"
#include "pointers.h"
#include "random_park.h"
//#include "STUBS/mpi.h"
#include "update.h"
#include "variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001
#define DELTA 1.005

// enum{PAIR,KSPACE,ATOM};
// enum{DIAMETER,CHARGE};

/* ---------------------------------------------------------------------- */

FixEPSExtract::FixEPSExtract(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec) error->all(FLERR,"Fix kinetics requires atom style bio");

  if (narg != 7) error->all(FLERR,"Illegal fix eps extract command: Missing arguments");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix eps extract command: calling steps should be positive integer");

  var = new char*[2];
  ivar = new int[2];

  int i;
  for (i = 0; i < 2; i++) {
    int n = strlen(&arg[4+i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i],&arg[4+i][2]);
  }

  seed = atoi(arg[6]);

  if (seed <= 0) error->all(FLERR,"Illegal fix eps extract command: seed should be greater than 0");

  // preExchangeCalled = false;

  // Random number generator, same for all procs
  random = new RanPark(lmp,seed);  

  find_maxid();

  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  // Set up renieghbouring here: required for re building the neighbour list: fix pour/ deposit
 
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  
}

/* ---------------------------------------------------------------------- */

FixEPSExtract::~FixEPSExtract()
{
  delete random;
  int i;
  for (i = 0; i < 2; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
}

/* ---------------------------------------------------------------------- */

int FixEPSExtract::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
   if need to restore per-atom quantities, create new fix STORE styles
------------------------------------------------------------------------- */


void FixEPSExtract::init()
{
     // fprintf(stdout, "called once?\n");
  if (!atom->radius_flag)
    error->all(FLERR,"Fix eps extract requires atom attribute diameter");

  int i;
  for (i = 0; i < 2; i++) {
    ivar[i] = input->variable->find(var[i]);
    if (ivar[i] < 0)
      error->all(FLERR,"Variable name for fix eps extract does not exist");
    if (!input->variable->equalstyle(ivar[i]))
      error->all(FLERR,"Variable for fix eps extract is invalid style");
  }

  int id = 0;
  for (int i = 1; i < group->ngroup; i++) {
    if (strcmp(group->names[i],"HET") == 0) {
      id = i;
      break;
    }
  }

  if (groupbit != pow(2, id)) {
    error->all(FLERR,"Fix eps extract is only valid for particles of type HET");
  }

  typeEPS = avec->typeEPS;

  if (typeEPS == 0) {
    error->all(FLERR,"Cannot find EPS type");
  }
}


void FixEPSExtract::post_integrate()
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  double EPSratio = input->variable->compute_equal(ivar[0]);
  double EPSdens = input->variable->compute_equal(ivar[1]);

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  for (i = 0; i < nall; i++) {
    if ((atom->mask[i] & groupbit) && atom->x[i][0] >= sublo[0] && atom->x[i][0] < subhi[0] &&
          atom->x[i][1] >= sublo[1] && atom->x[i][1] < subhi[1] &&
          atom->x[i][2] >= sublo[2] && atom->x[i][2] < subhi[2]) {
      // fprintf(stdout, "outerRadius/radius = %e\n", (outerRadius[i]/radius[i]));
      if ((avec->outerRadius[i]/atom->radius[i]) > EPSratio) {
        avec->outerMass[i] = (4.0*MY_PI/3.0)*(( avec->outerRadius[i] * avec->outerRadius[i]* avec->outerRadius[i])
      	    - (atom->radius[i]*atom->radius[i]*atom->radius[i]))*EPSdens;

        double splitF = 0.4 + (random->uniform()*0.2);

        double newOuterMass =  avec->outerMass[i] * splitF;
        double EPSMass = avec->outerMass[i] - newOuterMass;

        avec->outerMass[i] = newOuterMass;

        double density = atom->rmass[i] / (4.0*MY_PI/3.0 *
        								 atom->radius[i]*atom->radius[i]*atom->radius[i]);
        avec->outerRadius[i] = pow((3.0/(4.0*MY_PI))*((atom->rmass[i]/density)+( avec->outerMass[i]/EPSdens)),(1.0/3.0));

        double thetaD = random->uniform() * 2*MY_PI;
        double phiD = random->uniform() * (MY_PI);

        double oldX = atom->x[i][0];
        double oldY = atom->x[i][1];
        double oldZ = atom->x[i][2];

        //create child
        double childRadius = pow(((6*EPSMass)/(EPSdens*MY_PI)),(1.0/3.0))*0.5;
        double* coord = new double[3];
        double newX = oldX - ((childRadius+ avec->outerRadius[i])*cos(thetaD)*sin(phiD)*DELTA);
        double newY = oldY - ((childRadius+ avec->outerRadius[i])*sin(thetaD)*sin(phiD)*DELTA);
        double newZ = oldZ - ((childRadius+ avec->outerRadius[i])*cos(phiD)*DELTA);
        if (newX - childRadius < xlo) {
          newX = xlo + childRadius;
        }
        else if (newX + childRadius > xhi) {
          newX = xhi - childRadius;
        }
        if (newY - childRadius < ylo) {
          newY = ylo + childRadius;
        }
        else if (newY + childRadius > yhi) {
          newY = yhi - childRadius;
        }
        if (newZ - childRadius < zlo) {
          newZ = zlo + childRadius;
        }
        else if (newZ + childRadius > zhi) {
          newZ = zhi - childRadius;
        }
        coord[0] = newX;
        coord[1] = newY;
        coord[2] = newZ;
        find_maxid();
        atom->avec->create_atom(typeEPS,coord);
        // fprintf(stdout, "Created atom\n");
        int n = atom->nlocal - 1;
        atom->tag[n] = maxtag_all + 1;
        atom->mask[n] = avec->maskEPS;

        atom->v[n][0] = atom->v[i][0];
        atom->v[n][1] = atom->v[i][1];
        atom->v[n][2] = atom->v[i][2];
        atom->f[n][0] = atom->f[i][0];
        atom->f[n][1] = atom->f[i][1];
        atom->f[n][2] = atom->f[i][2];

        atom->omega[n][0] = atom->omega[i][0];
        atom->omega[n][1] = atom->omega[i][1];
        atom->omega[n][2] = atom->omega[i][2];
        atom->rmass[n] = EPSMass;
        avec->outerMass[n] = 0;

        atom->torque[n][0] = atom->torque[i][0];
        atom->torque[n][1] = atom->torque[i][1];
        atom->torque[n][2] = atom->torque[i][2];
        atom->radius[n] = childRadius;
        avec->outerRadius[n] = childRadius;
        //avec->atom_q[n] = 0;

        atom->natoms++;

        delete[] coord;
      }
    }
  }
  //fprintf(stdout, "Divided: %i\n", divided);

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;
}


/* ----------------------------------------------------------------------
   maxtag_all = current max atom ID for all atoms
------------------------------------------------------------------------- */


void FixEPSExtract::find_maxid()
{
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  maxtag_all = max;
  //MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
}




