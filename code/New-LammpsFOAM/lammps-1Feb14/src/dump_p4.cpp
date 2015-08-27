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
#include "dump_p4.h"
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

/* ---------------------------------------------------------------------- */

DumpP4::DumpP4(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump p4 command (5 arguments!)");
  if (binary || multiproc || multifile) error->all(FLERR,"DumpP4: this class does not accept 'binary', 'multiproc' or 'multifile' options");

  size_one = 10;
  size_conts = 8; // number of values per atom-atom contact
  size_walls = 7; // number of values per atom-wall contact
  sort_flag = 1;
  sortcol = 0;

  char *str = (char *) "%d  %d  %g  %g  %g  %g  %g  %g  %g  %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);

  str = (char *) "%d  %d  %g  %g  %g  %g  %g  %g\n"; // p1 p2 cx cy cz fx fy fz
  format4c = new char[strlen(str)];
  strcpy(format4c,str);

  str = (char *) "%d  %g  %g  %g  %g  %g  %g\n"; // p1 cx cy cz fx fy fz
  format4w = new char[strlen(str)];
  strcpy(format4w,str);

  maxbuf4c = 0;
  maxbuf4w = 0;
  buf4c = NULL;
  buf4w = NULL;
  nfix = 0;
}

/* ---------------------------------------------------------------------- */

DumpP4::~DumpP4()
{
  if(fp4p != NULL) fclose(fp4p);
  if(fp4c != NULL) fclose(fp4c);
  if(fp4w != NULL) fclose(fp4w);
}

/* ---------------------------------------------------------------------- */

void DumpP4::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // P4P file
  int len = strlen(filename) + 4;
  filename_p4p = new char[len];
  strcpy(filename_p4p,filename);
  strcat(filename_p4p,".p4p");

  // P4C file
  filename_p4c = new char[len];
  strcpy(filename_p4c,filename);
  strcat(filename_p4c,".p4c");

  int p4c = modify->find_compute("p4c");
  if (p4c<0)
    error->all(FLERR,"Illegal dump p4 call (compute p4c is required: [compute p4c all gral/local])");
  compute_pairs = modify->compute[p4c];


  nfix = 0;
  for (int i = 0; i < modify->nfix; i++)
  {
    if (modify->fix[i]->compute_local_flag)
      nfix++;
  }
  fix = new Fix*[nfix];

  nfix = 0;
  for (int i = 0; i < modify->nfix; i++)
  {
    if (modify->fix[i]->compute_local_flag)
    {
      fix[nfix] = modify->fix[i];
      nfix++;
    }
  }
  

  // P4W file
  filename_p4w = new char[len];
  strcpy(filename_p4w,filename);
  strcat(filename_p4w,".p4w");

  // open single file, one time only

  //if (multifile == 0) multifile always 0 with dump p4 style!
  openfile();
}

/* ---------------------------------------------------------------------- */

void DumpP4::write_header(bigint np)
{
  if (me == 0) {
    // P4P
    fprintf(fp4p,"\nTIMESTEP PARTICLES\n");
    fprintf(fp4p,BIGINT_FORMAT " " BIGINT_FORMAT "\n",update->ntimestep,np);
    fprintf(fp4p,"ID GROUP VOLUME MASS PX PY PZ VX VY VZ\n");
    // P4C
    fprintf(fp4c,"\nTIMESTEP CONTACTS\n");
    fprintf(fp4c,BIGINT_FORMAT " " BIGINT_FORMAT "\n",update->ntimestep,ctotal);
    fprintf(fp4c,"P1 P2 CX CY CZ FX FY FZ\n");
    // P4W
    fprintf(fp4w,"\nTIMESTEP CONTACTS\n");
    fprintf(fp4w,BIGINT_FORMAT " " BIGINT_FORMAT "\n",update->ntimestep,wtotal);
    fprintf(fp4w,"P1 CX CY CZ FX FY FZ\n");
  }
}

/* ---------------------------------------------------------------------- */

int DumpP4::count()
{
  //if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;

  double tmp = compute_pairs->compute_scalar(); // count contacts
  tmp = tmp;
  nconts = compute_pairs->size_local_rows;

  nwalls = 0;
  for (int i = 0; i < nfix; i++)
  {
    fix[i]->compute_local();
    nwalls += fix[i]->size_local_rows;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpP4::pack(int *ids)
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

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
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
      if (ids) ids[n++] = tag[i];
    }

  compute_pairs->compute_local();
  int nconts = compute_pairs->size_local_rows;

  m = 0;
  for (int i = 0; i < nconts; i++){
    p1 = (int) (compute_pairs->array_local[i][0]);
    p2 = (int) (compute_pairs->array_local[i][0]);
    dx = x[p1][0] - x[p2][0];
    dy = x[p1][1] - x[p2][1];
    dz = x[p1][2] - x[p2][2];
    dl = sqrt(dx*dx + dy*dy + dz*dz);
    dl = radius[p2] / dl;
    cx = x[p2][0] + dx * dl;
    cy = x[p2][1] + dy * dl;
    cz = x[p2][2] + dz * dl;

    buf4c[m++] = compute_pairs->array_local[i][0]; // tag1
    buf4c[m++] = compute_pairs->array_local[i][1]; // tag2
    buf4c[m++] = cx;
    buf4c[m++] = cy;
    buf4c[m++] = cz;
    buf4c[m++] = compute_pairs->array_local[i][2]; // Fx
    buf4c[m++] = compute_pairs->array_local[i][3]; // Fy
    buf4c[m++] = compute_pairs->array_local[i][4]; // Fz
    /*     for (int j = 0; j < size_conts; j++, m++) */
    /*       buf4c[m] = compute_pairs->array_local[i][j]; */
  }

  m = 0;
  for (int i = 0; i < nfix; i++)
  {
    fix[i]->compute_local();
    for (int j = 0; j < fix[i]->size_local_rows; j++){
      p1 = (int) fix[i]->array_local[j][0];

      buf4w[m++] = fix[i]->array_local[j][0]; // tag1
      
      buf4w[m++] = x[p1][0] + fix[i]->array_local[j][4]; // cx
      buf4w[m++] = x[p1][0] + fix[i]->array_local[j][5]; // cy
      buf4w[m++] = x[p1][0] + fix[i]->array_local[j][6]; // cz
      
      buf4w[m++] = fix[i]->array_local[j][1]; // fx
      buf4w[m++] = fix[i]->array_local[j][2]; // fy
      buf4w[m++] = fix[i]->array_local[j][3]; // fz
      
      /*       for (int k = 0; k < size_walls; k++, m++) */
      /*         buf4w[m] = fix[i]->array_local[j][k]; */
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpP4::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  if(multifile)
    error->one(FLERR,"DumpP4::openfile. Cannot use multifile option with P4 dump style");
    
  char *filecurrent_p4p;
  char *filecurrent_p4c;
  char *filecurrent_p4w;
  filecurrent_p4p = filename_p4p;
  filecurrent_p4c = filename_p4c;
  filecurrent_p4w = filename_p4w;

  // open one file on proc 0 or file on every proc

  if(multiproc)
    error->one(FLERR,"DumpP4::openfile. Cannot use multiproc option with P4 dump style");

  if (me == 0) {
    if (compressed) {
      error->one(FLERR,"DumpP4::openfile. Cannot use gzipped file with P4 dump style");
    } else if (binary) {
      error->one(FLERR,"DumpP4::openfile. Cannot use binary file with P4 dump style");
    } else if (append_flag) {
      fp4p = fopen(filecurrent_p4p,"a");
      fp4c = fopen(filecurrent_p4c,"a");
      fp4w = fopen(filecurrent_p4w,"a");
    } else {
      fp4p = fopen(filecurrent_p4p,"w");
      fp4c = fopen(filecurrent_p4c,"w");
      fp4w = fopen(filecurrent_p4w,"w");
    }

    if (fp4p == NULL) error->one(FLERR,"DumpP4::openfile. Cannot open P4P dump file");
    if (fp4c == NULL) error->one(FLERR,"DumpP4::openfile. Cannot open P4C dump file");
    if (fp4w == NULL) error->one(FLERR,"DumpP4::openfile. Cannot open P4W dump file");
  } else {
    fp4p = NULL;
    fp4c = NULL;
    fp4w = NULL;
  }

}

/* ---------------------------------------------------------------------- */

void DumpP4::write()
{
  // if file per timestep, open new file

  if (multifile) {
    //openfile();
    error->all(FLERR,"DumpP4::write. P4 dump file does not allow multifiles");
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
    error->all(FLERR,"DumpP4::write. P4 dump file does not allow multiprocessor files");
  else {
    MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
    
    MPI_Allreduce(&bnmec,&ctotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nmec,&nmaxc,1,MPI_INT,MPI_MAX,world);
    
    MPI_Allreduce(&bnmew,&wtotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nmew,&nmaxw,1,MPI_INT,MPI_MAX,world);
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
  if (nmaxc > maxbuf4c) {
    if ((bigint) nmaxc * size_conts > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf4c = nmaxc;
    memory->destroy(buf4c);
    memory->create(buf4c,maxbuf4c*size_conts,"dump:buf4c");
  }
  if (nmaxw > maxbuf4w) {
    if ((bigint) nmaxw * size_walls > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf4w = nmaxw;
    memory->destroy(buf4w);
    memory->create(buf4w,maxbuf4w*size_walls,"dump:buf4c");
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
    error->all(FLERR,"DumpP4::write. P4 dump file does not allow multiprocessor files");
  else {
    int tmp1,tmp2,tmp4,nlines,nlinec, nlinew;
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
          // p4c
          MPI_Irecv(buf4c,maxbuf4c*size_conts,MPI_DOUBLE,iproc,0,world,&request2);
          MPI_Send(&tmp2,0,MPI_INT,iproc,0,world);
          MPI_Wait(&request2,&status2);
          MPI_Get_count(&status2,MPI_DOUBLE,&nlinec);
          nlinec /= size_conts;
          // p4c
          MPI_Irecv(buf4w,maxbuf4w*size_walls,MPI_DOUBLE,iproc,0,world,&request3);
          MPI_Send(&tmp4,0,MPI_INT,iproc,0,world);
          MPI_Wait(&request3,&status3);
          MPI_Get_count(&status3,MPI_DOUBLE,&nlinew);
          nlinew /= size_walls;
        } else {
          nlines = nme;
          nlinec = nmec;
          nlinew = nmew;
        }

        write_data(nlines,buf);
        write_conts(nlinec,buf4c);
        write_walls(nlinew,buf4w);
      }
      if (flush_flag) {
        fflush(fp4p);
        fflush(fp4c);
        fflush(fp4w);
      }

    } else {
      MPI_Recv(&tmp1,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,world);
      
      MPI_Recv(&tmp2,0,MPI_INT,0,0,world,&status2);
      MPI_Rsend(buf4c,nmec*size_conts,MPI_DOUBLE,0,0,world);

      MPI_Recv(&tmp4,0,MPI_INT,0,0,world,&status3);
      MPI_Rsend(buf4w,nmew*size_walls,MPI_DOUBLE,0,0,world);
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpP4::write_data(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp4p,format,
        static_cast<int>(mybuf[m]),static_cast<int>(mybuf[m+1]),mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7],mybuf[m+8],mybuf[m+9]);
    m += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpP4::write_conts(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp4c,format4c,
        static_cast<int>(mybuf[m]),static_cast<int>(mybuf[m+1]),mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6],mybuf[m+7]);
    m += size_conts; // compute store [p1 p2 cx cy cz fx fy fz]
  }
}

/* ---------------------------------------------------------------------- */

void DumpP4::write_walls(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp4w,format4w,
        static_cast<int>(mybuf[m]),mybuf[m+1],mybuf[m+2],mybuf[m+3],mybuf[m+4],mybuf[m+5],mybuf[m+6]);
    m += size_walls; // fix store [p1 fx fy fz Dx Dy Dz]
  }
}

/* ---------------------------------------------------------------------- */
