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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_eps_extract.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "group.h"
#include "modify.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "kspace.h"
#include "fix_store.h"
#include "input.h"
#include "variable.h"
#include "respa.h"
#include "comm.h"
#include "domain.h"
#include "region.h"
#include "region_block.h"
#include "region_cylinder.h"
#include "random_park.h"
#include "math_extra.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001

// enum{PAIR,KSPACE,ATOM};
// enum{DIAMETER,CHARGE};

/* ---------------------------------------------------------------------- */

FixEPSExtract::FixEPSExtract(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
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
   

  // All gather arrays: fix pour
  // MPI_Comm_rank(world,&me);
  // MPI_Comm_size(world,&nprocs);
  // recvcounts = new int[nprocs];
  // displs = new int[nprocs];

  // Set up renieghbouring here: required for re building the neighbour list: fix pour/ deposit
 
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  



  // dynamic_group_allow = 1;
  // create_attribute = 1;

  // // count # of adaptations

  // ndiaadapt = 0;

  // int iarg = 4;
  // while (iarg < narg) {
  //   if (strcmp(arg[iarg],"pair") == 0) {
  //     if (iarg+6 > narg) error->all(FLERR,"Illegal fix diaadapt command");
  //     ndiaadapt++;
  //     iarg += 6;
  //   } else if (strcmp(arg[iarg],"kspace") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix diaadapt command");
  //     ndiaadapt++;
  //     iarg += 2;
  //   } else if (strcmp(arg[iarg],"atom") == 0) {
  //     if (iarg+3 > narg) error->all(FLERR,"Illegal fix diaadapt command");
  //     ndiaadapt++;
  //     iarg += 3;
  //   } else break;
  // }

  // if (nadapt == 0) error->all(FLERR,"Illegal fix adapt command");
  // adapt = new Adapt[nadapt];

  // // parse keywords

  // nadapt = 0;
  // diamflag = 0;
  // chgflag = 0;

  // iarg = 4;
  // while (iarg < narg) {
  //   if (strcmp(arg[iarg],"pair") == 0) {
  //     if (iarg+6 > narg) error->all(FLERR,"Illegal fix adapt command");
  //     adapt[nadapt].which = PAIR;
  //     int n = strlen(arg[iarg+1]) + 1;
  //     adapt[nadapt].pstyle = new char[n];
  //     strcpy(adapt[nadapt].pstyle,arg[iarg+1]);
  //     n = strlen(arg[iarg+2]) + 1;
  //     adapt[nadapt].pparam = new char[n];
  //     strcpy(adapt[nadapt].pparam,arg[iarg+2]);
  //     force->bounds(arg[iarg+3],atom->ntypes,
  //                   adapt[nadapt].ilo,adapt[nadapt].ihi);
  //     force->bounds(arg[iarg+4],atom->ntypes,
  //                   adapt[nadapt].jlo,adapt[nadapt].jhi);
  //     if (strstr(arg[iarg+5],"v_") == arg[iarg+5]) {
  //       n = strlen(&arg[iarg+5][2]) + 1;
  //       adapt[nadapt].var = new char[n];
  //       strcpy(adapt[nadapt].var,&arg[iarg+5][2]);
  //     } else error->all(FLERR,"Illegal fix adapt command");
  //     nadapt++;
  //     iarg += 6;
  //   } else if (strcmp(arg[iarg],"kspace") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
  //     adapt[nadapt].which = KSPACE;
  //     if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
  //       int n = strlen(&arg[iarg+1][2]) + 1;
  //       adapt[nadapt].var = new char[n];
  //       strcpy(adapt[nadapt].var,&arg[iarg+1][2]);
  //     } else error->all(FLERR,"Illegal fix adapt command");
  //     nadapt++;
  //     iarg += 2;
  //   } else if (strcmp(arg[iarg],"atom") == 0) {
  //     if (iarg+3 > narg) error->all(FLERR,"Illegal fix adapt command");
  //     adapt[nadapt].which = ATOM;
  //     if (strcmp(arg[iarg+1],"diameter") == 0) {
  //       adapt[nadapt].aparam = DIAMETER;
  //       diamflag = 1;
  //     } else if (strcmp(arg[iarg+1],"charge") == 0) {
  //       adapt[nadapt].aparam = CHARGE; 
  //       chgflag = 1; 
  //     } else error->all(FLERR,"Illegal fix adapt command");
  //     if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
  //       int n = strlen(&arg[iarg+2][2]) + 1;
  //       adapt[nadapt].var = new char[n];
  //       strcpy(adapt[nadapt].var,&arg[iarg+2][2]);
  //     } else error->all(FLERR,"Illegal fix adapt command");
  //     nadapt++;
  //     iarg += 3;
  //   } else break;
  // }

  // optional keywords

  // resetflag = 0;
  // scaleflag = 0;

  // while (iarg < narg) {
  //   if (strcmp(arg[iarg],"reset") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
  //     if (strcmp(arg[iarg+1],"no") == 0) resetflag = 0;
  //     else if (strcmp(arg[iarg+1],"yes") == 0) resetflag = 1;
  //     else error->all(FLERR,"Illegal fix adapt command");
  //     iarg += 2;
  //   } else if (strcmp(arg[iarg],"scale") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix adapt command");
  //     if (strcmp(arg[iarg+1],"no") == 0) scaleflag = 0;
  //     else if (strcmp(arg[iarg+1],"yes") == 0) scaleflag = 1;
  //     else error->all(FLERR,"Illegal fix adapt command");
  //     iarg += 2;
  //   } else error->all(FLERR,"Illegal fix adapt command");
  // }

  // allocate pair style arrays

  // int n = atom->ntypes;
  // for (int m = 0; m < nadapt; m++) {
  //   if (adapt[m].which == PAIR)
  //     memory->create(adapt[m].array_orig,n+1,n+1,"adapt:array_orig");
  // }

  // id_fix_diam = id_fix_chg = NULL;
}

/* ---------------------------------------------------------------------- */

FixEPSExtract::~FixEPSExtract()
{
  delete random;
  // delete [] recvcounts;
  // delete [] displs;
  int i;
  for (i = 0; i < 2; i++) {
    delete [] var[i];
  }
  delete [] var;
  delete [] ivar;
  // for (int m = 0; m < nadapt; m++) {
  //   delete [] adapt[m].var;
  //   if (adapt[m].which == PAIR) {
  //     delete [] adapt[m].pstyle;
  //     delete [] adapt[m].pparam;
  //     memory->destroy(adapt[m].array_orig);
  //   }
  // }
  // delete [] adapt;

  // // check nfix in case all fixes have already been deleted

  // if (id_fix_diam && modify->nfix) modify->delete_fix(id_fix_diam);
  // if (id_fix_chg && modify->nfix) modify->delete_fix(id_fix_chg);
  // delete [] id_fix_diam;
  // delete [] id_fix_chg;
}

/* ---------------------------------------------------------------------- */

int FixEPSExtract::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  // mask |= PRE_FORCE;
  // mask |= END_OF_STEP;
  return mask;
}

/* ----------------------------------------------------------------------
   if need to restore per-atom quantities, create new fix STORE styles
------------------------------------------------------------------------- */

// void FixAdapt::post_constructor()
// {
//   if (!resetflag) return;
//   if (!diamflag && !chgflag) return;

//   // new id = fix-ID + FIX_STORE_ATTRIBUTE
//   // new fix group = group for this fix

//   id_fix_diam = NULL;
//   id_fix_chg = NULL;

//   char **newarg = new char*[5];
//   newarg[1] = group->names[igroup];
//   newarg[2] = (char *) "STORE";
//   newarg[3] = (char *) "1";
//   newarg[4] = (char *) "1";

//   if (diamflag) {
//     int n = strlen(id) + strlen("_FIX_STORE_DIAM") + 1;
//     id_fix_diam = new char[n];
//     strcpy(id_fix_diam,id);
//     strcat(id_fix_diam,"_FIX_STORE_DIAM");
//     newarg[0] = id_fix_diam;
//     modify->add_fix(5,newarg);
//     fix_diam = (FixStore *) modify->fix[modify->nfix-1];

//     if (fix_diam->restart_reset) fix_diam->restart_reset = 0;
//     else {
//       double *vec = fix_diam->vstore;
//       double *radius = atom->radius;
//       int *mask = atom->mask;
//       int nlocal = atom->nlocal;

//       for (int i = 0; i < nlocal; i++) {
//         if (mask[i] & groupbit) vec[i] = radius[i];
//         else vec[i] = 0.0;
//       }
//     }
//   }

//   if (chgflag) {
//     int n = strlen(id) + strlen("_FIX_STORE_CHG") + 1;
//     id_fix_chg = new char[n];
//     strcpy(id_fix_chg,id);
//     strcat(id_fix_chg,"_FIX_STORE_CHG");
//     newarg[0] = id_fix_chg;
//     modify->add_fix(5,newarg);
//     fix_chg = (FixStore *) modify->fix[modify->nfix-1];

//     if (fix_chg->restart_reset) fix_chg->restart_reset = 0;
//     else {
//       double *vec = fix_chg->vstore;
//       double *q = atom->q;
//       int *mask = atom->mask;
//       int nlocal = atom->nlocal;

//       for (int i = 0; i < nlocal; i++) {
//         if (mask[i] & groupbit) vec[i] = q[i];
//         else vec[i] = 0.0;
//       }
//     }
//   }

//   delete [] newarg;
// }

/* ---------------------------------------------------------------------- */

// void FixDiaAdapt::end_of_step() {
// 	fprintf(stdout, "End called\n");
// 	preExchangeCalled = false;
// }

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
  if (groupbit != 2) {
    error->all(FLERR,"Fix eps extract is only valid for particles of type HET");
  }
}


// void FixDiaAdapt::pre_force(int vflag)
// {
//   if (nevery == 0) return;
//   if (update->ntimestep % nevery) return;
//  //  fprintf(stdout, "change_settings called\n");
//   change_dia();
//  //  fprintf(stdout, "before preexchange call\n");
//   // if (preExchangeCalled == false) {
//  	// pre_exchange();
//   // }
// }


void FixEPSExtract::pre_exchange()
{
  if (next_reneighbor != update->ntimestep) return;

  double EPSratio = input->variable->compute_equal(ivar[0]);
  double EPSdens = input->variable->compute_equal(ivar[1]);

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outerRadius = atom->outerRadius;
  double *outerMass = atom->outerMass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;
  // int nfix = modify->nfix;
  // Fix **fix = modify->fix;

  // int nnew = countNewAtoms(EPSratio);

  // find_maxid();

  double *sublo,*subhi;
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }




  // ncount = # of my atoms that overlap the insertion region
  // nprevious = total of ncount across all procs

  // int ncount = atom->nlocal;

  // int nprevious;
  // MPI_Allreduce(&ncount,&nprevious,1,MPI_INT,MPI_SUM,world);

  // // xmine is for my atoms
  // // xnear is for atoms from all procs + atoms to be inserted

  // double **xmine,**xnear;
  // memory->create(xmine,ncount,4,"fix_divide:xmine");
  // memory->create(xnear,nprevious+nnew,4,"fix_divide:xnear");
  // int nnear = nprevious;

  // // setup for allgatherv

  // int n = 4*ncount;
  // MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

  // displs[0] = 0;
  // for (int iproc = 1; iproc < nprocs; iproc++)
  //   displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  // // load up xmine array

  // double **x = atom->x;

  // ncount = 0;
  // for (i = 0; i < atom->nlocal; i++)
  // {
  //   xmine[ncount][0] = x[i][0];
  //   xmine[ncount][1] = x[i][1];
  //   xmine[ncount][2] = x[i][2];
  //   xmine[ncount][3] = radius[i];
  //   ncount++;
  // }

  // // perform allgatherv to acquire list of nearby particles on all procs

  // double *ptr = NULL;
  // if (ncount) ptr = xmine[0];
  // MPI_Allgatherv(ptr,4*ncount,MPI_DOUBLE,
  //                xnear[0],recvcounts,displs,MPI_DOUBLE,world);













  for (i = 0; i < nlocal; i++) {
    if ((mask[i] & groupbit) && atom->x[i][0] >= sublo[0] && atom->x[i][0] < subhi[0] &&
          atom->x[i][1] >= sublo[1] && atom->x[i][1] < subhi[1] &&
          atom->x[i][2] >= sublo[2] && atom->x[i][2] < subhi[2]) {
      // fprintf(stdout, "outerRadius/radius = %e\n", (outerRadius[i]/radius[i]));
      if ((outerRadius[i]/radius[i]) > EPSratio) {
        outerMass[i] = (4.0*MY_PI/3.0)*((outerRadius[i]*outerRadius[i]*outerRadius[i])-(radius[i]*radius[i]*radius[i]))*EPSdens;

        double splitF = 0.4 + (random->uniform()*0.2);

        double newOuterMass = outerMass[i] * splitF;
        double EPSMass = outerMass[i] - newOuterMass;

        outerMass[i] = newOuterMass;

        double density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
        outerRadius[i] = pow((3.0/(4.0*MY_PI))*((rmass[i]/density)+(outerMass[i]/EPSdens)),(1.0/3.0));

        double thetaD = random->uniform() * 2*MY_PI;
        double phiD = random->uniform() * (MY_PI);

        double oldX = atom->x[i][0];
        double oldY = atom->x[i][1];
        double oldZ = atom->x[i][2];


        //create child
        double childRadius = pow(((6*EPSMass)/(EPSdens*MY_PI)),(1.0/3.0))*0.5;
        double* coord = new double[3];
        double newX = oldX - (childRadius+outerRadius[i])*cos(thetaD)*sin(phiD);
        double newY = oldY - (childRadius+outerRadius[i])*sin(thetaD)*sin(phiD);
        double newZ = oldZ - (childRadius+outerRadius[i])*cos(phiD);
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
        atom->avec->create_atom(16,coord);
        // fprintf(stdout, "Created atom\n");
        int n = atom->nlocal - 1;
        atom->tag[n] = maxtag_all+1;
        atom->type[n] = 4;
        atom->mask[n] = 16;
        atom->v[n][0] = atom->v[i][0];
        atom->v[n][1] = atom->v[i][1];
        atom->v[n][2] = atom->v[i][2];
        atom->f[n][0] = atom->f[i][0];
        atom->f[n][1] = atom->f[i][1];
        atom->f[n][2] = atom->f[i][2];
        atom->omega[n][0] = atom->omega[i][0];
        atom->omega[n][1] = atom->omega[i][1];
        atom->omega[n][2] = atom->omega[i][2];
        rmass[n] = EPSMass;
        outerMass[n] = 0;
        atom->sub[n] = 0;
        atom->o2[n] = 0;
        atom->nh4[n] = 0;
        atom->no2[n] = 0;
        atom->no3[n] = 0;
        atom->torque[n][0] = atom->torque[i][0];
        atom->torque[n][1] = atom->torque[i][1];
        atom->torque[n][2] = atom->torque[i][2];
        radius[n] = childRadius;
        outerRadius[n] = childRadius;

        // for (j = 0; j < nfix; j++)
        //   if (fix[j]->create_attribute) fix[j]->set_arrays(n);

        // fprintf(stdout, "Set fields\n");
      //  modify->create_attribute(n);
       // fprintf(stdout, "Diameter of atom: %f\n", radius[n]*2);

        atom->natoms++;
      }
    }
  }

  //fprintf(stdout, "Divided: %i\n", divided);








    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }


  // free local memory

  // memory->destroy(xmine);
  // memory->destroy(xnear);

















  next_reneighbor += nevery;
}

// double FixDivide::getAverageMass()
// {
//   int *mask = atom->mask;
//   int nall = atom->nlocal + atom->nghost;
//   double *rmass = atom->rmass;
//   double averageMass = 0.0;
//   int numAtoms = 0;
//   int i;

//   for (i = 0; i < nall; i++) {
//     if (mask[i] & groupbit) {
//       averageMass += rmass[i];
//       numAtoms ++;
//     }
//   }

//   averageMass /= numAtoms;

//   return averageMass;
// }

// int FixEPSExtract::countNewAtoms(double EPSratio)
// {
//   int *mask = atom->mask;
//   int nall = atom->nlocal + atom->nghost;
//   double *radius = atom->radius;
//   double *outerRadius = atom->outerRadius;
//   int i;
//   int counter = 0;

//   for (i = 0; i < nall; i++) {
//     if (mask[i] & groupbit) {
//       if ((outerRadius[i]/radius[i]) > EPSratio) {
//         counter++;
//       }
//     }
//   }
//   return counter;
// }

/* ---------------------------------------------------------------------- */

// void FixDiaAdapt::pre_force_respa(int vflag, int ilevel, int)
// {
//   if (ilevel < nlevels_respa-1) return;
//   pre_force(vflag);
// }

/* ---------------------------------------------------------------------- */

// void FixDiaAdapt::post_run()
// {
//   if (resetflag) restore_settings();
// }


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
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
}


// int FixEPSExtract::overlap(int i)
// {
//   double delta;
//   delta = atom->radius[i] + radius_max;
//   double *boxlo = domain->boxlo;
//   double *boxhi = domain->boxhi;
//   double *prd = domain->prd;
//   int *periodicity = domain->periodicity;

//   double *x = atom->x[i];

//   if (domain->dimension == 3) {
//     if (region_style == 1) {
//       if (outside(0,x[0],xlo-delta,xhi+delta)) return 0;
//       if (outside(1,x[1],ylo-delta,yhi+delta)) return 0;
//       if (outside(2,x[2],lo_current-delta,hi_current+delta)) return 0;
//     } else {
//       double delx = x[0] - xc;
//       double dely = x[1] - yc;
//       double delz = 0.0;
//       domain->minimum_image(delx,dely,delz);
//       double rsq = delx*delx + dely*dely;
//       double r = rc + delta;
//       if (rsq > r*r) return 0;
//       if (outside(2,x[2],lo_current-delta,hi_current+delta)) return 0;
//     }
//   } else {
//     if (outside(0,x[0],xlo-delta,xhi+delta)) return 0;
//     if (outside(1,x[1],lo_current-delta,hi_current+delta)) return 0;
//   }

//   return 1;
// }


// int FixEPSExtract::outside(int dim, double value, double lo, double hi)
// {
//   double boxlo = domain->boxlo[dim];
//   double boxhi = domain->boxhi[dim];

//   if (domain->periodicity[dim]) {
//     if (lo < boxlo && hi > boxhi) {
//       return 0;
//     } else if (lo < boxlo) {
//       if (value > hi && value < lo + domain->prd[dim]) return 1;
//     } else if (hi > boxhi) {
//       if (value > hi - domain->prd[dim] && value < lo) return 1;
//     } else {
//       if (value < lo || value > hi) return 1;
//     }
//   } 

//   if (value < lo || value > hi) return 1;
//   return 0;
// }




