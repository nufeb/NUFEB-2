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

#include "string.h"
#include "dump_p3.h"
#include "atom.h"
#include "group.h"
#include "update.h"
#include "error.h"
#include "memory.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "compute.h"
#include "math_const.h" // MY_PI
#include "math.h"

using namespace LAMMPS_NS;

// customized by adding kaywords
enum{ATOM,COMPUTE,FIX};

enum{FX,FY,FZ,OMEGAX,OMEGAY,OMEGAZ,TQX,TQY,TQZ};

enum{P3,P4};

#define INVOKED_PERATOM 8


/* ---------------------------------------------------------------------- */
/*
 * 0: ID
 * 1: group-ID
 * 2: dump_type -> "p3"
 * 3: frequency
 * 4: outfilename
 * -- new [optionals]
 * 5: cont
 * 6: compute-ID (if "cont" is defined. "NULL" deactivate .p3c file)
 * 7: fix
 * 8: all | NULL | #nfix ("all" for all fix with compute_local_flag founded, or #nfix to defined how many fix will be defined, or "NULL" to deactivate .p3w file)
 * 9...: fix-IDs (if a number #nfix is defined previously)
 *10: vars
 *11: #vars
 *12: .....
 *13: style (output file style)
 *14: p3 | p4
*/
DumpP3::DumpP3(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump p3 command (5 arguments!)");
  if (binary || multiproc || multifile) error->all(FLERR,"DumpP3: this class does not accept 'binary', 'multiproc' or 'multifile' options");

  // store frequency of dump
  nevery = atoi(arg[3]);

  size_one = 10;
  size_conts = 5; // number of values per atom-atom contact
  size_walls = 7; // number of values per atom-wall contact
  sort_flag = 1;
  sortcol = 0;

  maxbuf3c = 0;
  maxbuf3w = 0;
  buf3c = NULL;
  buf3w = NULL;
  use_p3c = true;
  use_p3w = true;
  compute_pairs = NULL;
  //strcpy(compute_name,"p3c");
  nfix = 0;
  //fix_all = false;
  fix_all = true;
  nfix_names = 0;
  n_extra_vars = 0;
  header_extra = NULL;
  n_atom_vars = 0;
  var_atom = NULL;
  n_compute_vars = 0;
  var_compute = NULL;
  n_fix_vars = 0;
  var_fix = NULL;
  p3_style = P3;
  n_vcomputes = 0;
  n_vfixs = 0;

  bool cont_def = false;

  int iarg = 5;
  while (iarg<narg)
  {
    if ( strcmp(arg[iarg],"cont")==0 ) {
      if (iarg+1==narg)
        error->all(FLERR,"DumpP3: 'cont' keyword without argument in dump command");
      iarg++; // next argument
      if (strcmp(arg[iarg],"NULL")==0) {
        use_p3c = false;
      } else {
        compute_name = new char[strlen(arg[iarg])+1];
        strcpy(compute_name,arg[iarg]);
        cont_def = true;
      }
    } else if ( strcmp(arg[iarg],"fix")==0 ) {
      if (iarg+1==narg)
        error->all(FLERR,"DumpP3: 'fix' keyword without argument in dump command");

      iarg++;
      if (strcmp(arg[iarg],"all")==0) {
        fix_all = true;
      } else if (strcmp(arg[iarg],"NULL")==0) {
        use_p3w = false;
        fix_all = false;
      } else {
        nfix_names = atoi(arg[iarg]);
        fix_all = false;
      }
      if (nfix_names > 0){
        if (iarg+nfix_names>=narg)
          error->all(FLERR,"DumpP3: number of fix-ID defined higher than number of arguments in dump command");
        fix_name = new char*[nfix_names];
        for (int i=0; i<nfix_names; i++){
          iarg++;
          fix_name[i] = new char[strlen(arg[iarg])+1];
          strcpy(fix_name[i],arg[iarg]);
        }
      }

    } else if ( strcmp(arg[iarg],"vars")==0 ) {
      // definition of extra variables in .p3p file
      if (iarg+1==narg)
        error->all(FLERR,"DumpP3: 'var' keyword without argument in dump command");

      iarg++;
      if ( strcmp(arg[iarg],"NULL")==0) {
        n_extra_vars = 0;
      } else {
        n_extra_vars = atoi(arg[iarg]);

        if (n_extra_vars>0) {
          // read and parse extra per-atom variables
          if (iarg+n_extra_vars>=narg)
            error->all(FLERR,"DumpP3: number of extra variables higher that dump command arguments");
          iarg++;
          parse_extra_variables(n_extra_vars,arg+iarg);
          iarg+=n_extra_vars;
        } else {
          n_extra_vars = 0;
        }

      }

    } else if ( strcmp(arg[iarg],"style")==0 ) {
      // definition of P3 output file style (P3 or P4)
      if (iarg+1==narg)
        error->all(FLERR,"DumpP3: 'var' keyword without argument in dump command");
      iarg++;
      if ( strcmp(arg[iarg],"p3")==0 ) {
        p3_style = P3;
      } else if ( strcmp(arg[iarg],"p4")==0 ) {
        p3_style = P4;
      } else {
        error->all(FLERR,"DumpP3: 'style' argument not recognized in dump command");
      }

    } else {
      error->all(FLERR,"DumpP3: Invalid keyword in dump command");
    }
    iarg++;
  }

  if (use_p3c && !cont_def) {
    // use default name "p3c"
    compute_name = new char[4];
    strcpy(compute_name,"p3c");
  }

  // formats

  char *str = (char *) "%d  %d  %g  %g  %g  %g  %g  %g  %g  %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  str = (char *) "%d  %d  %g  %g  %g\n";
  format3c = new char[strlen(str)];
  strcpy(format3c,str);

  str = (char *) "%d  %g  %g  %g  %g  %g  %g\n";
  format3w = new char[strlen(str)];
  strcpy(format3w,str);

  // variables
  size_one += n_extra_vars; // add extra memory for variables
  if (p3_style==P3){
    size_conts = 5; // number of values per atom-atom contact
    size_walls = 7; // number of values per atom-wall contact
  } else { // P4
    size_conts = 8; // number of values per atom-atom contact
    size_walls = 7; // number of values per atom-wall contact
  }

}

/* ---------------------------------------------------------------------- */

DumpP3::~DumpP3()
{
  if(fp3p != NULL) fclose(fp3p);
  if(fp3c != NULL) fclose(fp3c);
  if(fp3w != NULL) fclose(fp3w);
}

/* ---------------------------------------------------------------------- */

void DumpP3::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + n_extra_vars*4 + 2;
  format = new char[n];
  strcpy(format,str);
  //for(int i=0; i<n_extra_vars; i++)
  //  strcat(format,"  %g");
  //strcat(format,"\n"); -> done in write_data()

  // P3P file
  int len = strlen(filename) + 4;
  filename_p3p = new char[len];
  strcpy(filename_p3p,filename);
  if (p3_style==P3)
    strcat(filename_p3p,".p3p");
  else
    strcat(filename_p3p,".p4p");

  // P3C file
  if (use_p3c)
  {
    //int p3c = modify->find_compute("p3c");
    int p3c = modify->find_compute(compute_name);
    if (p3c<0){
      //error->all(FLERR,"Illegal dump p3 call (compute p3c is required: [compute p3c all gral/local])");
      error->message(FLERR,"INFO: compute p3c (gran/local) not found. File .p3c/p4c not created");
      use_p3c = false;
    } else {
      use_p3c = true;
      compute_pairs = modify->compute[p3c];
      // prepare .p3c file
      filename_p3c = new char[len];
      strcpy(filename_p3c,filename);
      if (p3_style==P3)
        strcat(filename_p3c,".p3c");
      else
        strcat(filename_p3c,".p4c");
    }
  }

  // P3W file
  nfix = 0;
  if (fix_all){
    for (int i = 0; i < modify->nfix; i++) {
      if (modify->fix[i]->compute_local_flag)
        nfix++;
    }
  } else if (use_p3w) {
    int ifix = 0;
    for (int i = 0; i < nfix_names; i++) {
      ifix = modify->find_fix(fix_name[i]);
      if (ifix>0 && modify->fix[ifix]->compute_local_flag)
        nfix++;
      if (ifix<0)
        error->message(FLERR,"INFO: DumpP3: fix-ID not found.");
    }
  }

  if (nfix>0)
    fix = new Fix*[nfix];

  nfix = 0;
  if (fix_all) {
    for (int i = 0; i < modify->nfix; i++){
      if (modify->fix[i]->compute_local_flag){
        fix[nfix] = modify->fix[i];
        nfix++;
      }
    }
  } else if (use_p3w) {
    int ifix = 0;
    for (int i = 0; i < nfix_names; i++) {
      ifix = modify->find_fix(fix_name[i]);
      if (ifix>0 && modify->fix[ifix]->compute_local_flag){
        fix[nfix] = modify->fix[ifix];
        nfix++;
      }
    }
  }

  if (nfix>0)
  {
    use_p3w = true;
    filename_p3w = new char[len];
    strcpy(filename_p3w,filename);
    if (p3_style==P3)
      strcat(filename_p3w,".p3w");
    else
      strcat(filename_p3w,".p4w");
  } else if (use_p3w) {
    use_p3w = false;
    error->message(FLERR,"INFO: fix (wall/gran) not found. File .p3w/p4w not created");
  }

  // open single file, one time only

  //if (multifile == 0) multifile always 0 with dump p3 style!
  openfile();
}

/* ---------------------------------------------------------------------- */

void DumpP3::write_header(bigint np)
{
  if (me == 0) {
    // P3P
    fprintf(fp3p,"\nTIMESTEP PARTICLES\n");
    fprintf(fp3p,BIGINT_FORMAT " " BIGINT_FORMAT "\n",update->ntimestep,np);
    //fprintf(fp3p,"ID GROUP RADIUS MASS PX PY PZ VX VY VZ\n");
    if (p3_style==P3)
      fprintf(fp3p,"ID GROUP RADIUS MASS PX PY PZ VX VY VZ");
    else
      fprintf(fp3p,"ID GROUP VOLUME MASS PX PY PZ VX VY VZ");
    if (header_extra)
      fprintf(fp3p,"%s",header_extra);
    fprintf(fp3p,"\n");
    // P3C
    if (use_p3c)
    {
      fprintf(fp3c,"\nTIMESTEP CONTACTS\n");
      fprintf(fp3c,BIGINT_FORMAT " " BIGINT_FORMAT "\n",update->ntimestep,ctotal);
      if (p3_style==P3)
        fprintf(fp3c,"P1 P2 FX FY FZ\n");
      else // P4
        fprintf(fp3c,"P1 P2 CX CY CZ FX FY FZ\n");
    }
    // P3W
    if (use_p3w)
    {
      fprintf(fp3w,"\nTIMESTEP CONTACTS\n");
      fprintf(fp3w,BIGINT_FORMAT " " BIGINT_FORMAT "\n",update->ntimestep,wtotal);
      if (p3_style==P3)
        fprintf(fp3w,"P1 FX FY FZ NX NY NZ\n");
      else
        fprintf(fp3w,"P1 CX CY CZ FX FY FZ\n");
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpP3::split_var_name(char* label, char* var, char* alias, int& type)
{
  type = ATOM;

  int len = strlen(label);
  if ( strncmp(label,"c_",2)==0 ) {

    type = COMPUTE;
    char *ptr = strchr(label,':'); // verify if exist alias
    if (ptr){
      strcpy(alias,ptr+1); // label = c_ID[#]:alias
      // for(len=2;len<strlen(label);len++) if (label[len]==':') break;
      len = (ptr-label);  // need to be verified!!! signed and unsigned int comparison...
    } else
      strcpy(alias,label+2);

    strncpy(var,label,len);
    var[len]='\0'; // force string end
    //len -= 2;

  } else if ( strncmp(label,"f_",2)==0 ) {
    type = FIX;

    char *ptr = strchr(label,':'); // verify if exist alias
    if (ptr){
      strcpy(alias,ptr+1); // label = var:alias
      len = (ptr-label);  // need to be verified!!! signed and unsigned int comparison...
    } else
      strcpy(alias,label+2);

    strncpy(var,label,len);
    var[len]='\0'; // force string end
  
  } else { // atom
    type = ATOM;

    char *ptr = strchr(label,':'); // verify if exist alias
    if (ptr){
      strcpy(alias,ptr+1); // label = var:alias
      len = (ptr-label);  // need to be verified!!! signed and unsigned int comparison...
    } else
      strcpy(alias,label);

    strncpy(var,label,len);
    var[len]='\0'; // force string end

  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::add_compute(char* var, int& i_compute, int& i_var)
{

  // process c_var[]
  char* ptr = strchr(var,'[');
  int len = strlen(var);
  if (ptr) {
    if (var[len-1] != ']')
      error->all(FLERR,"DumpP3: Invalid compute attribute in dump custom command");
    var_index[i_var] = atoi(ptr+1);
    *ptr = '\0';
  } else {
    var_index[i_var] = 0;
  }

  int id;
  for(id=0; id<n_vcomputes; id++)
    if (strcmp(id_vcomputes[id],var+2)==0) break;
  if (id<n_vcomputes) {
    var_compute[i_compute] = id;
  } else {

    id_vcomputes = (char **) memory->srealloc(id_vcomputes,(n_vcomputes+1)*sizeof(char *),"DumpP3:id_vcomputes");
    
    //delete [] vcomputes;
    //vcomputes = new Compute*[n_vcomputes+1];

    int n = strlen(var+2) + 1;
    id_vcomputes[n_vcomputes] = new char[n];
    strcpy(id_vcomputes[n_vcomputes],var+2);
    n_vcomputes++;

    var_compute[i_compute] = n_vcomputes-1;
  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::add_fix(char* var, int& i_fix, int& i_var)
{

  // process f_var[]
  char* ptr = strchr(var,'[');
  int len = strlen(var);
  if (ptr) {
    if (var[len-1] != ']')
      error->all(FLERR,"DumpP3: Invalid fix attribute in dump custom command");
    var_index[i_var] = atoi(ptr+1);
    *ptr = '\0';
  } else {
    var_index[i_var] = 0;
  }

  int id;
  for(id=0; id<n_vfixs; id++)
    if (strcmp(id_vfixs[id],var+2)==0) break;
  if (id<n_vfixs) {
    var_fix[i_fix] = id;
  } else {

    id_vfixs = (char **) memory->srealloc(id_vfixs,(n_vfixs+1)*sizeof(char *),"DumpP3:id_vfixs");
    
    int n = strlen(var+2) + 1;
    id_vfixs[n_vfixs] = new char[n];
    strcpy(id_vfixs[n_vfixs],var+2);
    n_vfixs++;

    var_fix[i_fix] = n_vfixs-1;
  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::parse_extra_variables(int narg,char **arg)
{

  char var[64];
  char name[64];
  int  type;

  var_type = new int[narg];
  var_ptr  = new double*[narg];
  var_stride = new int[narg];
  var_index = new int[narg];

  // count variables
  for(int i=0; i<narg; i++)
  {
    split_var_name(arg[i],var,name,type);
    if ( type==COMPUTE ) {
      n_compute_vars++;
    } else if ( type==FIX ) {
      n_fix_vars++;
    } else if ( type==ATOM ) {
      n_atom_vars++;
    } else {
      error->all(FLERR,"DumpP3: extra variable in dump command not recognized");
    }
    var_type[i] = type;
  }
  
  // fill variables

  if (n_atom_vars)
    var_atom = new int[n_atom_vars];
  if (n_compute_vars)
    var_compute = new int[n_compute_vars];
  if (n_fix_vars)
    var_fix = new int[n_fix_vars];

  char buffer[1024] = "\0"; // temporal
  char *string;
  string = buffer;

  int n_atoms = 0;
  int n_computes = 0;
  int n_fixs = 0;
  for(int i=0; i<narg; i++)
  {
    split_var_name(arg[i],var,name,type);
    if ( strncmp(var,"c_",2)==0 ) {
      add_compute(var,n_computes,i);
      n_computes++;
    } else if ( strncmp(var,"f_",2)==0 ) {
      add_fix(var,n_fixs,i);
      n_fixs++;
    } else if ( strcmp(var,"omegax")==0 ) {
      if (!atom->omega_flag)
        error->all(FLERR,"DumpP3: Atom property not available: omega");
      var_atom[n_atoms++] = OMEGAX;
    } else if (strcmp(var,"omegay")==0 ) {
      if (!atom->omega_flag)
        error->all(FLERR,"DumpP3: Atom property not available: omega");
      var_atom[n_atoms++] = OMEGAY;
    } else if (strcmp(var,"omegaz")==0 ) {
      if (!atom->omega_flag)
        error->all(FLERR,"DumpP3: Atom property not available: omega");
      var_atom[n_atoms++] = OMEGAZ;
    } else if (strcmp(var,"fx")==0 ) {
      var_atom[n_atoms++] = FX;
    } else if (strcmp(var,"fy")==0 ) {
      var_atom[n_atoms++] = FY;
    } else if (strcmp(var,"fz")==0 ) {
      var_atom[n_atoms++] = FZ;
    } else if (strcmp(var,"tqx")==0 ) {
      if (!atom->torque_flag)
        error->all(FLERR,"DumpP3: Atom property not available: torque");
      var_atom[n_atoms++] = TQX;
    } else if (strcmp(var,"tqy")==0 ) {
      if (!atom->torque_flag)
        error->all(FLERR,"DumpP3: Atom property not available: torque");
      var_atom[n_atoms++] = TQY;
    } else if (strcmp(var,"tqz")==0 ) {
      if (!atom->torque_flag)
        error->all(FLERR,"DumpP3: Atom property not available: torque");
      var_atom[n_atoms++] = TQZ;
    } else {
      error->all(FLERR,"DumpP3: extra variable in dump command not recognized");
    }
    strcat(string," ");
    strcat(string,name);
  }

  header_extra = new char[strlen(string)+2];
  strcpy(header_extra,string);

  // Allocate computes
  if (n_vcomputes){
    vcomputes = new Compute*[n_vcomputes];
    for(int i=0;i<n_vcomputes;i++){
      int id = modify->find_compute(id_vcomputes[i]);
      if (id<0)
        error->all(FLERR,"DumpP3: compute-ID not found");
      vcomputes[i] = modify->compute[id];
    }
  }

  // Allocate fixes
  if (n_vfixs){
    vfixs = new Fix*[n_vfixs];
    for(int i=0;i<n_vfixs;i++){
      int id = modify->find_fix(id_vfixs[i]);
      if (id<0) 
        error->all(FLERR,"DumpP3: fix-ID not found");
      if (nevery % modify->fix[id]->peratom_freq)
        error->all(FLERR,"DumpP3: dump and fix not computed at compatible times");
      vfixs[i] = modify->fix[id];
    }

    int i_fix = 0;
    for(int i=0; i<n_extra_vars; i++){
      int id = var_fix[i_fix];
      if (var_type[i]==FIX){
        // verify column
        if (vfixs[id]->peratom_flag == 0)
          error->all(FLERR,"DumpP3: fix does not compute per-atom info");
        if (var_index[i] == 0 && vfixs[id]->size_peratom_cols > 0)
          error->all(FLERR,"DumpP3: fix does not compute per-atom vector");
        if (var_index[i] > 0 && vfixs[id]->size_peratom_cols == 0)
          error->all(FLERR,"DumpP3: fix does not compute per-atom array");
        if (var_index[i] > 0 &&  var_index[i] > vfixs[id]->size_peratom_cols)
          error->all(FLERR,"DumpP3: fix vector is accessed out-of-range");
        i_fix++;
      }
    }

  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::init_extra_variables_pointer()
{

  int n_atoms = 0;
  int i_compute = 0;
  int i_fixs = 0;
  int icomp, ifix;

  // invoke Computes for per-atom quantities
  if (n_vcomputes) {
    for (int i = 0; i < n_vcomputes; i++)
      if (!(vcomputes[i]->invoked_flag & INVOKED_PERATOM)) {
        vcomputes[i]->compute_peratom();
        vcomputes[i]->invoked_flag |= INVOKED_PERATOM;
      }
  }

  for (int i=0; i<n_extra_vars; i++)
  {
    if (var_type[i]==COMPUTE) {
     
      icomp = var_compute[i_compute];
      if (var_index[i] == 0) {
        var_ptr[i] = vcomputes[icomp]->vector_atom;
        var_stride[i] = 1;
      } else {
        var_ptr[i] = &vcomputes[icomp]->array_atom[0][var_index[i]-1];
        var_stride[i] = vcomputes[icomp]->size_peratom_cols;
      }

      i_compute++;

    } else if (var_type[i]==FIX) {

      ifix = var_fix[i_fixs];
      if (var_index[i] == 0) {
        var_ptr[i] = vfixs[ifix]->vector_atom;
        var_stride[i] = 1;
      } else {
        var_ptr[i] = &vfixs[ifix]->array_atom[0][var_index[i]-1];
        var_stride[i] = vfixs[ifix]->size_peratom_cols;
      }
      i_fixs++;

    } else if (var_type[i]==ATOM) {

      if (var_atom[n_atoms]==OMEGAX) {
        var_ptr[i]    = &atom->omega[0][0];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==OMEGAY) {
        var_ptr[i]    = &atom->omega[0][1];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==OMEGAZ) {
        var_ptr[i]    = &atom->omega[0][2];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==FX) {
        var_ptr[i]    = &atom->f[0][0];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==FY) {
        var_ptr[i]    = &atom->f[0][1];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==FZ) {
        var_ptr[i]    = &atom->f[0][2];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==TQX) {
        var_ptr[i]    = &atom->torque[0][0];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==TQY) {
        var_ptr[i]    = &atom->torque[0][1];
        var_stride[i] = 3;
      } else if (var_atom[n_atoms]==TQZ) {
        var_ptr[i]    = &atom->torque[0][2];
        var_stride[i] = 3;

      }
      n_atoms++;

    }

  }

}

/* ---------------------------------------------------------------------- */

int DumpP3::count()
{
  //if (igroup == 0) return atom->nlocal;
  double dummy;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;

  if (use_p3c)
  {
    dummy = compute_pairs->compute_scalar(); // count contacts
    nconts = compute_pairs->size_local_rows;
  } else
    nconts = 0;

  nwalls = 0;
  for (int i = 0; i < nfix; i++) // nfix > 0 if use_p3w
  {
    fix[i]->compute_local();
    nwalls += fix[i]->size_local_rows;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpP3::pack(int *ids)
{
  if (p3_style==P3)
    pack_style_p3(ids);
  else // P4
    pack_style_p4(ids);
}

/* ---------------------------------------------------------------------- */

void DumpP3::pack_style_p3(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->rmass;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  // particles
  init_extra_variables_pointer();

  m = n = 0;
  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = radius[i];
      buf[m++] = mass[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];
      for (int j=0; j<n_extra_vars; j++)
        buf[m++] = *var_ptr[j];
      if (ids) ids[n++] = tag[i];
    }
    for (int j=0; j<n_extra_vars; j++)
      var_ptr[j] += var_stride[j];
  }

  // particle contacts

  int nconts = 0;
  if (use_p3c)
  {
    compute_pairs->compute_local();
    nconts = compute_pairs->size_local_rows;
  }

  m = 0;
  for (int i = 0; i < nconts; i++)
    for (int j = 0; j < size_conts; j++)
      buf3c[m++] = compute_pairs->array_local[i][j];

  // wall contacts

  m = 0;
  for (int i = 0; i < nfix; i++){
    fix[i]->compute_local();
    for (int j = 0; j < fix[i]->size_local_rows; j++)
      for (int k = 0; k < size_walls; k++) 
        buf3w[m++] = fix[i]->array_local[j][k];
  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::pack_style_p4(int *ids)
{
  int m,n, p1, p2;
  double cx, cy, cz, dx, dy, dz, dl, vol;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->rmass;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  static const double cte = MathConst::MY_PI * 4.0 / 3.0;

  // particles
  init_extra_variables_pointer();

  m = n = 0;
  for (int i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) {
      vol = cte * radius[i]*radius[i]*radius[i];
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = vol; //radius[i];
      buf[m++] = mass[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];
      for (int j=0; j<n_extra_vars; j++)
        buf[m++] = *var_ptr[j];
      if (ids) ids[n++] = tag[i];
    }
    for (int j=0; j<n_extra_vars; j++)
      var_ptr[j] += var_stride[j];
  }

  // particle contacts

  int nconts = 0;
  if (use_p3c)
  {
    compute_pairs->compute_local();
    nconts = compute_pairs->size_local_rows;
  }

  m = 0;
  for (int i = 0; i < nconts; i++){
    p1 = (int) (compute_pairs->array_local[i][0]);
    p2 = (int) (compute_pairs->array_local[i][1]);
    dx = x[p1][0] - x[p2][0];
    dy = x[p1][1] - x[p2][1];
    dz = x[p1][2] - x[p2][2];
    dl = sqrt(dx*dx + dy*dy + dz*dz);
    dl = radius[p2] / dl;
    cx = x[p2][0] + dx * dl;
    cy = x[p2][1] + dy * dl;
    cz = x[p2][2] + dz * dl;

    buf3c[m++] = compute_pairs->array_local[i][0]; // tag1
    buf3c[m++] = compute_pairs->array_local[i][1]; // tag2
    buf3c[m++] = cx;
    buf3c[m++] = cy;
    buf3c[m++] = cz;
    buf3c[m++] = compute_pairs->array_local[i][2]; // Fx
    buf3c[m++] = compute_pairs->array_local[i][3]; // Fy
    buf3c[m++] = compute_pairs->array_local[i][4]; // Fz
  }

  // wall contacts

  m = 0;
  for (int i = 0; i < nfix; i++)
  {
    fix[i]->compute_local();
    for (int j = 0; j < fix[i]->size_local_rows; j++){
      p1 = (int) fix[i]->array_local[j][0];

      buf3w[m++] = fix[i]->array_local[j][0]; // tag1
      
      buf3w[m++] = x[p1][0] + fix[i]->array_local[j][4]; // cx
      buf3w[m++] = x[p1][0] + fix[i]->array_local[j][5]; // cy
      buf3w[m++] = x[p1][0] + fix[i]->array_local[j][6]; // cz
      
      buf3w[m++] = fix[i]->array_local[j][1]; // fx
      buf3w[m++] = fix[i]->array_local[j][2]; // fy
      buf3w[m++] = fix[i]->array_local[j][3]; // fz
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  if(multifile)
    error->one(FLERR,"DumpP3::openfile. Cannot use multifile option with P3 dump style");
    
  char *filecurrent_p3p;
  char *filecurrent_p3c;
  char *filecurrent_p3w;
  filecurrent_p3p = filename_p3p;
  filecurrent_p3c = filename_p3c;
  filecurrent_p3w = filename_p3w;

  // open one file on proc 0 or file on every proc

  if(multiproc)
    error->one(FLERR,"DumpP3::openfile. Cannot use multiproc option with P3 dump style");

  if (me == 0) {
    if (compressed) {
      error->one(FLERR,"DumpP3::openfile. Cannot use gzipped file with P3 dump style");
    } else if (binary) {
      error->one(FLERR,"DumpP3::openfile. Cannot use binary file with P3 dump style");
    } else if (append_flag) {
      fp3p = fopen(filecurrent_p3p,"a");
      if(use_p3c) fp3c = fopen(filecurrent_p3c,"a");
      if(use_p3w) fp3w = fopen(filecurrent_p3w,"a");
    } else {
      fp3p = fopen(filecurrent_p3p,"w");
      if(use_p3c) fp3c = fopen(filecurrent_p3c,"w");
      if(use_p3w) fp3w = fopen(filecurrent_p3w,"w");
    }

    if (fp3p == NULL) error->one(FLERR,"DumpP3::openfile. Cannot open P3P dump file");
    if (use_p3c && fp3c == NULL) error->one(FLERR,"DumpP3::openfile. Cannot open P3C dump file");
    if (use_p3w && fp3w == NULL) error->one(FLERR,"DumpP3::openfile. Cannot open P3W dump file");
  } else {
    fp3p = NULL;
    fp3c = NULL;
    fp3w = NULL;
  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::write()
{
  // if file per timestep, open new file

  if (multifile) {
    //openfile();
    error->all(FLERR,"DumpP3::write. P3 dump file does not allow multifiles");
  }

  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }

  // nme = # of dump lines this proc will contribute to dump

  nme = count();
  bigint bnme = nme;

  int nmec = nconts;
  bigint bnmec = nmec;

  int nmew = nwalls;
  bigint bnmew = nmew;

  // ntotal = total # of dump lines
  // nmax = max # of dump lines on any proc

  int nmax, nmaxc, nmaxw;
  if (multiproc)
    //nmax = nme;
    error->all(FLERR,"DumpP3::write. P3 dump file does not allow multiprocessor files");
  else {
    MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
   
    if (use_p3c) {
      MPI_Allreduce(&bnmec,&ctotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
      MPI_Allreduce(&nmec,&nmaxc,1,MPI_INT,MPI_MAX,world);
    }
    
    if (use_p3w) {
      MPI_Allreduce(&bnmew,&wtotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
      MPI_Allreduce(&nmew,&nmaxw,1,MPI_INT,MPI_MAX,world);
    }
  }


  // write timestep header

  if (!multiproc)
    write_header(ntotal);

  // insure proc 0 can receive everyone's info
  // limit nmax*size_one to int since used as arg in MPI_Rsend() below
  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }
  if (use_p3c && nmaxc > maxbuf3c) {
    if ((bigint) nmaxc * size_conts > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf3c = nmaxc;
    memory->destroy(buf3c);
    memory->create(buf3c,maxbuf3c*size_conts,"dump:buf3c");
  }
  if (use_p3w && nmaxw > maxbuf3w) {
    if ((bigint) nmaxw * size_walls > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf3w = nmaxw;
    memory->destroy(buf3w);
    memory->create(buf3w,maxbuf3w*size_walls,"dump:buf3c");
  }
  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  if (sort_flag && sortcol == 0)
    pack(ids);
  else
    pack(NULL);
  if (sort_flag)
    sort();

  // multiproc = 1 = each proc writes own data to own file 
  // multiproc = 0 = all procs write to one file thru proc 0
  //   proc 0 pings each proc, receives it's data, writes to file
  //   all other procs wait for ping, send their data to proc 0

  if (multiproc)
    //write_data(nme,buf);
    error->all(FLERR,"DumpP3::write. P3 dump file does not allow multiprocessor files");
  else {
    int tmp1,tmp2,tmp3,nlines,nlinec, nlinew;
    MPI_Status status, status2, status3;
    MPI_Request request, request2, request3;

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
        if (iproc) {
          MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,iproc,0,world,&request);
          MPI_Send(&tmp1,0,MPI_INT,iproc,0,world);
          MPI_Wait(&request,&status);
          MPI_Get_count(&status,MPI_DOUBLE,&nlines);
          nlines /= size_one;
          // p3c
          if (use_p3c) {
            MPI_Irecv(buf3c,maxbuf3c*size_conts,MPI_DOUBLE,iproc,0,world,&request2);
            MPI_Send(&tmp2,0,MPI_INT,iproc,0,world);
            MPI_Wait(&request2,&status2);
            MPI_Get_count(&status2,MPI_DOUBLE,&nlinec);
            nlinec /= size_conts;
          }
          // p3c
          if (use_p3w) {
            MPI_Irecv(buf3w,maxbuf3w*size_walls,MPI_DOUBLE,iproc,0,world,&request3);
            MPI_Send(&tmp3,0,MPI_INT,iproc,0,world);
            MPI_Wait(&request3,&status3);
            MPI_Get_count(&status3,MPI_DOUBLE,&nlinew);
            nlinew /= size_walls;
          }
        } else {
          nlines = nme;
          nlinec = nmec;
          nlinew = nmew;
        }

        write_data(nlines,buf);
        if(use_p3c) write_conts(nlinec,buf3c);
        if(use_p3w) write_walls(nlinew,buf3w);
      }
      if (flush_flag) {
        fflush(fp3p);
        if(use_p3c) fflush(fp3c);
        if(use_p3w) fflush(fp3w);
      }

    } else {
      MPI_Recv(&tmp1,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,world);
     
      if (use_p3c) {
        MPI_Recv(&tmp2,0,MPI_INT,0,0,world,&status2);
        MPI_Rsend(buf3c,nmec*size_conts,MPI_DOUBLE,0,0,world);
      }

      if (use_p3w) {
        MPI_Recv(&tmp3,0,MPI_INT,0,0,world,&status3);
        MPI_Rsend(buf3w,nmew*size_walls,MPI_DOUBLE,0,0,world);
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpP3::write_data(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp3p,format,
        static_cast<int>(mybuf[m]),static_cast<int>(mybuf[m+1]),mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8],mybuf[m+9]);
    for(int j = 10; j< size_one; j++)
      fprintf(fp3p,"  %g",mybuf[m+j]);
    fprintf(fp3p,"\n");
    m += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpP3::write_conts(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp3c,format3c,
        static_cast<int>(mybuf[m]),static_cast<int>(mybuf[m+1]),mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_conts; // compute store [p1 p2 fx fy fz energ]
  }
}

/* ---------------------------------------------------------------------- */

void DumpP3::write_walls(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp3w,format3w,
        static_cast<int>(mybuf[m]),mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6]);
    m += size_walls; // fix store [p1 fx fy fz Dx Dy Dz]
  }
}

/* ---------------------------------------------------------------------- */
