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
#include "fix_diaadapt.h"
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
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

// enum{PAIR,KSPACE,ATOM};
// enum{DIAMETER,CHARGE};

/* ---------------------------------------------------------------------- */

FixDiaAdapt::FixDiaAdapt(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix diaadapt command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix diaadapt command");

  int n = strlen(&arg[4][2]) + 1;
  var = new char[n];
  strcpy(var,&arg[4][2]);

  growthFactor = atof(arg[5]);

  preExchangeCalled = false;

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

FixDiaAdapt::~FixDiaAdapt()
{
  delete [] var;
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

int FixDiaAdapt::setmask()
{
  int mask = 0;
  //mask |= PRE_EXCHANGE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
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

void FixDiaAdapt::end_of_step() {
	fprintf(stdout, "End called\n");
	preExchangeCalled = false;
}

void FixDiaAdapt::init()
{
  if (!atom->radius_flag)
    error->all(FLERR,"Fix diaadapt requires atom attribute diameter");

  ivar = input->variable->find(var);
  if (ivar < 0)
    error->all(FLERR,"Variable name for fix adapt does not exist");
  if (!input->variable->equalstyle(ivar))
    error->all(FLERR,"Variable for fix adapt is invalid style");




 //  int i,j;

 //  // allow a dynamic group only if ATOM attribute not used

 //  if (group->dynamic[igroup])
 //    for (int i = 0; i < nadapt; i++)
 //      if (adapt[i].which == ATOM) 
 //        error->all(FLERR,"Cannot use dynamic group with fix adapt atom");
  
 //  // setup and error checks

 //  anypair = 0;

 //  for (int m = 0; m < nadapt; m++) {
 //    Adapt *ad = &adapt[m];

 //    ad->ivar = input->variable->find(ad->var);
 //    if (ad->ivar < 0)
 //      error->all(FLERR,"Variable name for fix adapt does not exist");
 //    if (!input->variable->equalstyle(ad->ivar))
 //      error->all(FLERR,"Variable for fix adapt is invalid style");

 //    if (ad->which == PAIR) {
 //      anypair = 1;
 //      Pair *pair = NULL;

 //      if (lmp->suffix_enable) {
 //        char psuffix[128];
 //        strcpy(psuffix,ad->pstyle);
 //        strcat(psuffix,"/");
 //        strcat(psuffix,lmp->suffix);
 //        pair = force->pair_match(psuffix,1);
 //      }
 //      if (pair == NULL) pair = force->pair_match(ad->pstyle,1);
 //      if (pair == NULL) error->all(FLERR,"Fix adapt pair style does not exist");
 //      void *ptr = pair->extract(ad->pparam,ad->pdim);
 //      if (ptr == NULL) 
 //        error->all(FLERR,"Fix adapt pair style param not supported");

 //      ad->pdim = 2;
 //      if (ad->pdim == 0) ad->scalar = (double *) ptr;
 //      if (ad->pdim == 2) ad->array = (double **) ptr;

 //      // if pair hybrid, test that ilo,ihi,jlo,jhi are valid for sub-style

 //      if (ad->pdim == 2 && (strcmp(force->pair_style,"hybrid") == 0 ||
 //                            strcmp(force->pair_style,"hybrid/overlay") == 0)) {
 //        PairHybrid *pair = (PairHybrid *) force->pair;
 //        for (i = ad->ilo; i <= ad->ihi; i++)
 //          for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
 //            if (!pair->check_ijtype(i,j,ad->pstyle))
 //              error->all(FLERR,"Fix adapt type pair range is not valid for "
 //                         "pair hybrid sub-style");
 //      }

 //    } else if (ad->which == KSPACE) {
 //      if (force->kspace == NULL)
 //        error->all(FLERR,"Fix adapt kspace style does not exist");
 //      kspace_scale = (double *) force->kspace->extract("scale");

 //    } else if (ad->which == ATOM) {
 //      if (ad->aparam == DIAMETER) {
 //        if (!atom->radius_flag)
 //          error->all(FLERR,"Fix adapt requires atom attribute diameter");
 //        int nall = atom->nlocal + atom->nghost;
 //        double *radius = atom->radius;
 //        memory->create(adapt[m].radius_orig,nall,"adapt:radius_orig");
 //        for (i = 0; i < nall; i++) {
 //          ad->radius_orig[i] = radius[i];
 //        }
 //      }
 //      if (ad->aparam == CHARGE) {
	// if (!atom->q_flag)
	//   error->all(FLERR,"Fix adapt requires atom attribute charge");
 //      }
 //    }
 //  }

 //  // make copy of original pair array values

 //  for (int m = 0; m < nadapt; m++) {
 //    Adapt *ad = &adapt[m];
 //    if (ad->which == PAIR && ad->pdim == 2) {
 //      for (i = ad->ilo; i <= ad->ihi; i++)
 //        for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
 //          ad->array_orig[i][j] = ad->array[i][j];
 //    }
 //  }

 //  // fixes that store initial per-atom values
  
 //  if (id_fix_diam) {
 //    int ifix = modify->find_fix(id_fix_diam);
 //    if (ifix < 0) error->all(FLERR,"Could not find fix adapt storage fix ID");
 //    fix_diam = (FixStore *) modify->fix[ifix];
 //  }
 //  if (id_fix_chg) {
 //    int ifix = modify->find_fix(id_fix_chg);
 //    if (ifix < 0) error->all(FLERR,"Could not find fix adapt storage fix ID");
 //    fix_chg = (FixStore *) modify->fix[ifix];
 //  }

 //  if (strstr(update->integrate_style,"respa"))
 //    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

// void FixDiaAdapt::setup_pre_force(int vflag)
// {
//   change_settings();
// }

/* ---------------------------------------------------------------------- */

// void FixDiaAdapt::setup_pre_force_respa(int vflag, int ilevel)
// {
//   if (ilevel < nlevels_respa-1) return;
//   setup_pre_force(vflag);
// }

/* ---------------------------------------------------------------------- */

void FixDiaAdapt::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  // fprintf(stdout, "change_settings called\n");
  change_settings();
  // fprintf(stdout, "In between\n");
  if (preExchangeCalled == false) {
  	pre_exchange();
  }
}


void FixDiaAdapt::pre_exchange()
{
  preExchangeCalled = true;
  double density;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;
  double averageMass = 0.0;
  int numAtoms = 0;


  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      averageMass += rmass[i];
      numAtoms ++;
    }
  }

  averageMass /= numAtoms;

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      if (rmass[i] >= growthFactor*averageMass) {
        double splitF = 0.5;//0.3 + (random->uniform() * 0.4);
        double parentMass = rmass[i] * splitF;
        double childMass = rmass[i] - parentMass;

        double thetaD = MY_PI/2;//random->uniform() * MY_PI;
        double phiD = MY_PI/4;//random->uniform() * (MY_PI/2);

        double oldX = atom->x[i][0];
        double oldY = atom->x[i][1];
        double oldZ = atom->x[i][2];


        //Update parent
        rmass[i] = parentMass;
        radius[i] = pow(((6*rmass[i])/(density*MY_PI)),(1.0/3.0))*0.5;
        atom->x[i][0] = oldX + radius[i]*cos(thetaD)*sin(phiD);
        atom->x[i][1] = oldY + radius[i]*sin(thetaD)*sin(phiD);
        atom->x[i][2] = oldZ + radius[i]*cos(phiD);
        fprintf(stdout, "Diameter of atom: %f\n", radius[i]*2);

        // fprintf(stdout, "Moved and resized parent\n");

        //create child
        double childRadius = pow(((6*childMass)/(density*MY_PI)),(1.0/3.0))*0.5;
        double* coord = new double[3];
        coord[0] = oldX - childRadius*cos(thetaD)*sin(phiD);
        coord[1] = oldY - childRadius*sin(thetaD)*sin(phiD);
        coord[2] = oldZ - childRadius*cos(phiD);
        atom->avec->create_atom(mask[i],coord);
        // fprintf(stdout, "Created atom\n");
        int n = atom->nlocal - 1;
        atom->tag[n] = n+1;
        atom->type[n] = atom->type[i];
        atom->mask[n] = mask[i];
        atom->image[n] = atom->image[i];
        atom->v[n][0] = atom->v[i][0];
        atom->v[n][1] = atom->v[i][1];
        atom->v[n][2] = atom->v[i][2];
        rmass[n] = childMass;
        radius[n] = childRadius;
        // fprintf(stdout, "Set fields\n");
        modify->create_attribute(n);
        fprintf(stdout, "Diameter of atom: %f\n", radius[n]*2);

        atom->natoms++;
      }
    }
  }

  fprintf(stdout, "Number of atoms: %i\n", atom->natoms);
  fprintf(stdout, "pre_exchange called\n");
}

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
   change pair,kspace,atom parameters based on variable evaluation
------------------------------------------------------------------------- */

void FixDiaAdapt::change_settings()
{
  // int i,j;

  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // for (int m = 0; m < nadapt; m++) {
  //   Adapt *ad = &adapt[m];
  double value = input->variable->compute_equal(ivar);

    // set global scalar or type pair array values

    // if (ad->which == PAIR) {
    //   if (ad->pdim == 0) {
    //     if (scaleflag) *ad->scalar = value * ad->scalar_orig;
    //     else *ad->scalar = value;
    //   } else if (ad->pdim == 2) {
    //     if (scaleflag)
    //       for (i = ad->ilo; i <= ad->ihi; i++)
    //         for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
    //           ad->array[i][j] = value*ad->array_orig[i][j];
    //     else
    //       for (i = ad->ilo; i <= ad->ihi; i++)
    //         for (j = MAX(ad->jlo,i); j <= ad->jhi; j++)
    //           ad->array[i][j] = value;
    //   }

    // // set kspace scale factor

    // } else if (ad->which == KSPACE) {
    //   *kspace_scale = value;

    // // set per atom values, also make changes for ghost atoms

    // } else if (ad->which == ATOM) {

      // reset radius from diameter
      // also scale rmass to new value

      // if (ad->aparam == DIAMETER) {
    // int mflag = 0;
    // if (atom->rmass_flag) mflag = 1;
  double density;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int i;

  for (i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      density = rmass[i] / (4.0*MY_PI/3.0 *
                      radius[i]*radius[i]*radius[i]);
      double oldMass = rmass[i];
      rmass[i] = rmass[i]*(1 + (value*nevery));
      if (rmass[i] <= 0) {
        rmass[i] = oldMass;
      }
      radius[i] = pow(((6*rmass[i])/(density*MY_PI)),(1.0/3.0))*0.5;
    }
  }

    // if (mflag == 0) {
    //   for (i = 0; i < nall; i++)
    //     if (mask[i] & groupbit)
    //       radius[i] = 0.5*value;
    // } else {
    //   for (i = 0; i < nall; i++)
    //     if (mask[i] & groupbit) {
    //       density = rmass[i] / (4.0*MY_PI/3.0 *
    //                           radius[i]*radius[i]*radius[i]);
    //       radius[i] = 0.5*value;
    //       rmass[i] = 4.0*MY_PI/3.0 *
    //           radius[i]*radius[i]*radius[i] * density;
    //     }
    // }

  modify->addstep_compute(update->ntimestep + nevery);


  // re-initialize pair styles if any PAIR settings were changed
  // this resets other coeffs that may depend on changed values,
  // and also offset and tail corrections

  // if (anypair) force->pair->reinit();

  // // reset KSpace charges if charges have changed

  // if (chgflag && force->kspace) force->kspace->qsum_qsq();
}

/* ----------------------------------------------------------------------
   restore pair,kspace,atom parameters to original values
------------------------------------------------------------------------- */

// void FixAdapt::restore_settings()
// {
//   for (int m = 0; m < nadapt; m++) {
//     Adapt *ad = &adapt[m];
//     if (ad->which == PAIR) {
//       if (ad->pdim == 0) *ad->scalar = ad->scalar_orig;
//       else if (ad->pdim == 2) {
//         for (int i = ad->ilo; i <= ad->ihi; i++)
//           for (int j = MAX(ad->jlo,i); j <= ad->jhi; j++)
//             ad->array[i][j] = ad->array_orig[i][j];
//       }

//     } else if (ad->which == KSPACE) {
//       *kspace_scale = 1.0;

//     } else if (ad->which == ATOM) {
//       if (diamflag) {
//         double density;

//         double *vec = fix_diam->vstore;
//         double *radius = atom->radius;
//         double *rmass = atom->rmass;
//         int *mask = atom->mask;
//         int nlocal = atom->nlocal;

//         for (int i = 0; i < nlocal; i++)
//           if (mask[i] & groupbit) {
//             density = rmass[i] / (4.0*MY_PI/3.0 *
//                                   radius[i]*radius[i]*radius[i]);
//             radius[i] = vec[i];
//             rmass[i] = 4.0*MY_PI/3.0 * radius[i]*radius[i]*radius[i] * density;
//           }
//       }
//       if (chgflag) {
//         double *vec = fix_chg->vstore;
//         double *q = atom->q;
//         int *mask = atom->mask;
//         int nlocal = atom->nlocal;

//         for (int i = 0; i < nlocal; i++)
//           if (mask[i] & groupbit) q[i] = vec[i];
//       }
//     }
//   }

//   if (anypair) force->pair->reinit();
//   if (chgflag && force->kspace) force->kspace->qsum_qsq();
// }

// /* ----------------------------------------------------------------------
//    initialize one atom's storage values, called when atom is created
// ------------------------------------------------------------------------- */

// void FixAdapt::set_arrays(int i)
// {
//   if (fix_diam) fix_diam->vstore[i] = atom->radius[i];
//   if (fix_chg) fix_chg->vstore[i] = atom->q[i];
// }