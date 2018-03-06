/*
 * fix_metabolism.h
 *
 *  Created on: 15 Aug 2016
 *      Author: bowen
 */

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

#include "fix_bio_kinetics.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm> 

#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"

#include "bio.h"
#include "atom_vec_bio.h"
#include "fix_bio_kinetics_ph.h"
#include "fix_bio_kinetics_thermo.h"
#include "pointers.h"
#include "variable.h"
#include "modify.h"
#include "update.h"
#include "force.h"
#include "group.h"
#include "domain.h"
#include "fix_bio_kinetics_diffusion.h"
#include "fix_bio_kinetics_energy.h"
#include "fix_bio_kinetics_monod.h"
#include "comm.h"
#include "fix_bio_fluid.h"

#ifdef OUTPUT_GRID
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

using namespace std;

#define BUFMIN 1000

/* ---------------------------------------------------------------------- */

FixKinetics::FixKinetics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");

  bio = avec->bio;

  if (narg < 15)
    error->all(FLERR, "Not enough arguments in fix kinetics command");

  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix kinetics command: calling steps should be positive integer");

  if (strcmp(arg[15],"no") == 0) demflag = 0;
  else if (strcmp(arg[15],"yes") == 0) demflag = 1;
  else error->all(FLERR,"Illegal demflag parameter (yes or no)");

  var = new char*[7];
  ivar = new int[7];

  nx = atoi(arg[4]);
  ny = atoi(arg[5]);
  nz = atoi(arg[6]);

//  nout = 0;
//  if (narg > 14) { // grid output
//    nout = force->inumeric(FLERR, arg[14]);
//  }

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  } else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  for (int i = 0; i < 7; i++) {
    int n = strlen(&arg[7 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[7 + i][2]);
  }

  for (int i = 0; i < 3; i++) {
    // considering that the grid will always have a cubic cell (i.e. stepx = stepy = stepz)
    subnlo[i] = domain->sublo[i] / stepz;
    subnhi[i] = domain->subhi[i] / stepz;
    sublo[i] = subnlo[i] * stepz;
    subhi[i] = subnhi[i] * stepz;
    subn[i] = subnhi[i] - subnlo[i];
  }

  bnz = domain->boxhi[2] / stepz + 1;
  maxheight = domain->boxhi[2];
}

/* ---------------------------------------------------------------------- */

FixKinetics::~FixKinetics() {
  int i;
  for (i = 0; i < 7; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(gYield);
  memory->destroy(activity);
  memory->destroy(nuR);
  memory->destroy(nuS);
  memory->destroy(nuBS);
  memory->destroy(qGas);
  memory->destroy(DRGCat);
  memory->destroy(DRGAn);
  memory->destroy(kEq);
  memory->destroy(Sh);

  memory->destroy(recvcells);
  memory->destroy(sendcells);
  memory->destroy(recvbegin);
  memory->destroy(recvend);
  memory->destroy(sendbegin);
  memory->destroy(sendend);
  memory->destroy(fV);

  delete[] nuConv;
}

/* ---------------------------------------------------------------------- */

int FixKinetics::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init() {
  for (int n = 0; n < 7; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix kinetics does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix kinetics is invalid style");
  }

  // register fix kinetics with this class
  diffusion = NULL;
  energy = NULL;
  ph = NULL;
  thermo = NULL;
  monod = NULL;
  nufebFoam = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics/growth/energy") == 0) {
      energy = static_cast<FixKineticsEnergy *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/diffusion") == 0) {
      diffusion = static_cast<FixKineticsDiffusion *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/ph") == 0) {
      ph = static_cast<FixKineticsPH *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/thermo") == 0) {
      thermo = static_cast<FixKineticsThermo *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "kinetics/growth/monod") == 0) {
      monod = static_cast<FixKineticsMonod *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "nufebFoam") == 0) {
      nufebFoam = static_cast<FixFluid *>(lmp->modify->fix[j]);
    }
  }

  if (bio->nnus == 0)
    error->all(FLERR, "fix_kinetics requires # of Nutrients inputs");
  else if (bio->nuGCoeff == NULL && energy != NULL)
    error->all(FLERR, "fix_kinetics requires Nutrient Energy inputs");
  else if (bio->iniS == NULL)
    error->all(FLERR, "fix_kinetics requires Nutrients inputs");

  ngrids = subn[0] * subn[1] * subn[2];

  nnus = bio->nnus;
  int ntypes = atom->ntypes;

  temp = input->variable->compute_equal(ivar[0]);
  rth = input->variable->compute_equal(ivar[1]);
  gVol = input->variable->compute_equal(ivar[2]);
  gasTrans = input->variable->compute_equal(ivar[3]);
  iph = input->variable->compute_equal(ivar[4]);
  diffT = input->variable->compute_equal(ivar[5]);
  bl = input->variable->compute_equal(ivar[6]);

  nuConv = new int[nnus + 1]();
  nuS = memory->create(nuS, nnus + 1, ngrids, "kinetics:nuS");
  nuR = memory->create(nuR, nnus + 1, ngrids, "kinetics:nuR");
  nuBS = memory->create(nuBS, nnus + 1, "kinetics:nuBS");
  qGas = memory->create(qGas, nnus + 1, ngrids, "kinetics:nuGas");
  gYield = memory->create(gYield, ntypes + 1, ngrids, "kinetic:gYield");
  activity = memory->create(activity, nnus + 1, 5, ngrids, "kinetics:activity");
  DRGCat = memory->create(DRGCat, ntypes + 1, ngrids, "kinetics:DRGCat");
  DRGAn = memory->create(DRGAn, ntypes + 1, ngrids, "kinetics:DRGAn");
  kEq = memory->create(kEq, nnus + 1, 4, "kinetics:kEq");
  Sh = memory->create(Sh, ngrids, "kinetics:Sh");
  fV = memory->create(fV, 3, ngrids, "kinetcis:fV");

  //initialize grid yield, inlet concentration, consumption
  for (int j = 0; j < ngrids; j++) {
    fV[0][j] = 0;
    fV[1][j] = 0;
    fV[2][j] = 0;
    for (int i = 1; i <= ntypes; i++) {
      gYield[i][j] = bio->yield[i];
      DRGCat[i][j] = 0;
      DRGAn[i][j] = 0;
    }
    for (int i = 1; i <= nnus; i++) {
      nuS[i][j] = bio->iniS[i][0];
      nuR[i][j] = 0;
      qGas[i][j] = 0;
    }
  }

  reset_isConv();
  update_bgrids();

  if (energy != NULL) {
    init_keq();
    init_activity();
  }

  recv_buff_size = BUFMIN;
  send_buff_size = BUFMIN;
  recvcells = memory->create(recvcells, recv_buff_size, "kinetics::recvcells");
  sendcells = memory->create(sendcells, send_buff_size, "kinetics::sendcells");
  recvbegin = memory->create(recvbegin, comm->nprocs, "kinetics::recvbegin");
  recvend = memory->create(recvend, comm->nprocs, "kinetics::recvend");
  sendbegin = memory->create(sendbegin, comm->nprocs, "kinetics::sendbegin");
  sendend = memory->create(sendend, comm->nprocs, "kinetics::sendend");

  for (int i = 0; i < comm->nprocs; i++)
  {
    recvbegin[i] = 0;
    recvend[i] = 0;
    sendbegin[i] = 0;
    sendend[i] = 0;
  }

  // Fitting initial domain decomposition to the grid 
  for (int i = 0; i < comm->procgrid[0]; i++) {
    int n = nx * i * 1.0 / comm->procgrid[0];
    comm->xsplit[i] = (double)n / nx;
  }
  for (int i = 0; i < comm->procgrid[1]; i++) {
    int n = ny * i * 1.0 / comm->procgrid[1];
    comm->ysplit[i] = (double)n / ny;
  }
  for (int i = 0; i < comm->procgrid[2]; i++) {
    int n = nz * i * 1.0 / comm->procgrid[2];
    comm->zsplit[i] = (double)n / nz;
  }
  domain->set_local_box();

  borders();
}

/* ---------------------------------------------------------------------- */

void FixKinetics::borders()
{
  if (comm->nprocs < 2 || diffusion == NULL) return;

  // communicate grid extent
  int *gridlo = new int[3 * comm->nprocs];
  MPI_Allgather(&subnlo[0], 3, MPI_INT, gridlo, 3, MPI_INT, world);
  int *gridhi = new int[3 * comm->nprocs];
  MPI_Allgather(&subnhi[0], 3, MPI_INT, gridhi, 3, MPI_INT, world);

  int nrecv = 0;
  int nsend = 0;
  // look for intersections
  Grid grid(subnlo, subnhi);
  Grid extgrid = extend(grid);
  for (int p = 0; p < comm->nprocs; p++) {
    recvbegin[p] = nrecv;
    sendbegin[p] = nsend;
    Grid other(&gridlo[3 * p], &gridhi[3 * p]);
    if (p != comm->me) {
      // identify which cell we are going to receive
      send_recv_cells(extgrid, grid, other, nsend, nrecv);
      // check for periodic boundary conditions
      // TODO: check if diffusion is NULL
      if (extgrid.lower[0] < 0 && diffusion->xbcflag == 0) {
	send_recv_cells(extgrid, grid, translate(other, -nx, 0, 0), nsend, nrecv);
      }
      if (extgrid.upper[0] > nx && diffusion->xbcflag == 0) {      
	send_recv_cells(extgrid, grid, translate(other, nx, 0, 0), nsend, nrecv);
      }
      if (extgrid.lower[1] < 0 && diffusion->ybcflag == 0) {
	send_recv_cells(extgrid, grid, translate(other, 0, -ny, 0), nsend, nrecv);
      }
      if (extgrid.upper[1] > ny && diffusion->ybcflag == 0) {      
	send_recv_cells(extgrid, grid, translate(other, 0, ny, 0), nsend, nrecv);
      }
      if (extgrid.lower[2] < 0 && diffusion->zbcflag == 0) {
	send_recv_cells(extgrid, grid, translate(other, 0, 0, -nz), nsend, nrecv);
      }
      if (extgrid.upper[2] > nz && diffusion->zbcflag == 0) {      
	send_recv_cells(extgrid, grid, translate(other, 0, 0, nz), nsend, nrecv);
      }
    }
    recvend[p] = nrecv;
    sendend[p] = nsend;
  }

  delete [] gridlo;
  delete [] gridhi;
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init_keq() {
  // water Kj/mol
  double dG0H2O = -237.18;
  int nnus = bio->nnus;
  double **nuGCoeff = bio->nuGCoeff;

  for (int i = 1; i < nnus + 1; i++) {
    kEq[i][0] = exp((dG0H2O + nuGCoeff[i][0] - nuGCoeff[i][1]) / (-rth * temp));
    for (int j = 1; j < 4; j++) {
      double coeff = 0.0;

      if (nuGCoeff[i][j + 1] > 10000) {
        coeff = j * 10001;
      } else {
        coeff = 0;
      }

      kEq[i][j] = exp((nuGCoeff[i][j + 1] + coeff - nuGCoeff[i][j]) / (-rth * temp));
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixKinetics::init_activity() {
  int nnus = bio->nnus;
  double *denm = memory->create(denm, nnus + 1, "kinetics:denm");
  double gSh = pow(10, -iph);

  for (int k = 1; k < nnus + 1; k++) {
    for (int j = 0; j < ngrids; j++) {
      double iniNuS = bio->iniS[k][0];
      Sh[j] = gSh;
      denm[k] = (1 + kEq[k][0]) * gSh * gSh * gSh + kEq[k][1] * gSh * gSh + kEq[k][2] * kEq[k][3] * gSh + kEq[k][3] * kEq[k][2] * kEq[k][1];
      if (denm[k] == 0) {
        lmp->error->all(FLERR, "denm returns a zero value");
      }
      // not hydrated form acitivity
      activity[k][0][j] = kEq[k][0] * nuS[k][j] * gSh * gSh * gSh / denm[k];
      // fully protonated form activity
      activity[k][1][j] = nuS[k][j] * gSh * gSh * gSh / denm[k];
      // 1st deprotonated form activity
      activity[k][2][j] = nuS[k][j] * gSh * gSh * kEq[k][1] / denm[k];
      // 2nd deprotonated form activity
      activity[k][3][j] = nuS[k][j] * gSh * kEq[k][1] * kEq[k][2] / denm[k];
      // 3rd deprotonated form activity
      activity[k][4][j] = nuS[k][j] * kEq[k][1] * kEq[k][2] * kEq[k][3] / denm[k];

      if (strcmp(bio->nuName[k], "h") == 0) {
        activity[k][1][j] = gSh;
      }
    }
  }

  memory->destroy(denm);
}

/* ---------------------------------------------------------------------- */

void FixKinetics::pre_force(int vflag) {
  bool flag = true;
  if (nufebFoam != NULL)
    demflag = demflag || nufebFoam->demflag;

  if (nevery == 0)
    flag = false;
  if (update->ntimestep % nevery)
    flag = false;
  if(demflag)
    flag = false;

  if (flag)
    integration();

#ifdef OUTPUT_GRID
  if (update->ntimestep % nout == 0) {
    for (int i = 1; i <= nnus; i++)
      output_nutrient_info(i);
    if (energy != NULL) {
      for (int i = 1; i <= atom->ntypes; i++)
	output_bacteria_info(i);
    }
  }
#endif
}

/* ----------------------------------------------------------------------
   main kenetics integration loop
 ------------------------------------------------------------------------- */

void FixKinetics::integration() {
  int iteration = 0;
  bool isConv = false;

  gflag = 0;
  update_bgrids();
  reset_nuR();

  while (!isConv) {
    iteration++;
    isConv = true;

    // solve for reaction term, no growth happens here
    if (energy != NULL) {
      if (ph != NULL) ph->solve_ph();
      else init_activity();

      thermo->thermo();
      energy->growth(diffT, gflag);
    } else if (monod != NULL) {
      monod->growth(diffT, gflag);
    }

    // solve for diffusion and advection
    if (diffusion != NULL) {
      nuConv = diffusion->diffusion(nuConv, iteration, diffT);
    } else {
      reset_nuR();
      break;
    }

    // check for convergence
    for (int i = 1; i <= nnus; i++) {
      if (!nuConv[i]) {
        isConv = false;
        break;
      }
    }

    if (iteration >= 10000) isConv = true;
  }

  if (comm->me == 0) printf( "number of iteration: %i \n", iteration);

  gflag = 1;
  reset_isConv();

  // microbe growth
  if (energy != NULL) energy->growth(update->dt*nevery, gflag);
  if (monod != NULL) monod->growth(update->dt*nevery, gflag);

  // solve mass balance of nutrients in bulk liquid
  // the concentration of o2 in the bulk liquid is kept constant by aeration
  if (diffusion != NULL && diffusion->bulkflag == 1) {
    for (int i = 1; i <= nnus; i++) {
      if (strcmp(bio->nuName[i], "o2") != 0) {
        diffusion->compute_bulk(i);
      }
    }
  }

  if (diffusion != NULL) diffusion->update_nuS();
}

/* ----------------------------------------------------------------------
   get maximum height of biofilm
 ------------------------------------------------------------------------- */

double FixKinetics::getMaxHeight() {
  const int nlocal = atom->nlocal;
  double * const * const x = atom->x;
  double * const r = atom->radius;
  double maxh = 0;

  for (int i=0; i < nlocal; i++) {
    if((x[i][2] + r[i]) > maxh)
      maxh = x[i][2] + r[i]*2;
  }

  double global_max;
  MPI_Allreduce(&maxh, &global_max, 1, MPI_DOUBLE, MPI_MAX, world);

  return global_max;
}

void FixKinetics::update_bgrids() {
  if (bl > 0) {
    maxheight = getMaxHeight();
    bnz = (int)((bl + maxheight) / stepz) + 1;
    bgrids = subn[0] * subn[1] * MIN(subn[2], MAX(0, bnz - subnlo[2]));
  }
  else {
    bgrids = subn[0] * subn[1] * subn[2];
  }
}

bool FixKinetics::is_inside(int i) {
  if (atom->x[i][0] < sublo[0] || atom->x[i][0] >= subhi[0] ||
      atom->x[i][1] < sublo[1] || atom->x[i][1] >= subhi[1] ||
      atom->x[i][2] < sublo[2] || atom->x[i][2] >= subhi[2])
    return false;
  return true;
}  

int FixKinetics::position(int i) {
  // get index of grid containing i
  int xpos = (atom->x[i][0] - sublo[0]) / stepz;
  int ypos = (atom->x[i][1] - sublo[1]) / stepz;
  int zpos = (atom->x[i][2] - sublo[2]) / stepz;
  int pos = xpos + ypos * subn[0] + zpos * subn[0] * subn[1];

  if (pos >= bgrids) {
    printf("Too big! pos=%d   size = %i\n", pos, bgrids);
  }

  return pos;
}

void FixKinetics::add_cells(const Grid &basegrid, const Grid &grid, int *cells, int index)
{
  for (int k = grid.lower[2]; k < grid.upper[2]; k++) {
    for (int j = grid.lower[1]; j < grid.upper[1]; j++) {
      for (int i = grid.lower[0]; i < grid.upper[0]; i++) {
        cells[index++] = get_linear_index(basegrid, i, j, k);
      }
    }
  }
}

bool FixKinetics::is_intersection_valid(const Grid &g)
{
  // check if the intersection is empty
  if (is_empty(g))
    return false;
  // check if the intersection is a corner
  int n[3];
  for (int i = 0; i < 3; i++)
    n[i] = g.upper[i] - g.lower[i];
  if ((n[0] > 1 && n[1] > 1) || (n[0] > 1 && n[2] > 1) || (n[1] > 1 && n[2] > 1))
    return true;
  return false;
}

void FixKinetics::send_recv_cells(const Grid &basegrid, const Grid &grid, const Grid &other, int &nsend, int &nrecv)
{
  // identify which cells we need to recv
  Grid recvgrid = intersect(extend(grid), other);
  int n = cell_count(recvgrid);
  if (recv_buff_size < nrecv + n) { // not enough memory to fit cells
    recv_buff_size += (n / BUFMIN + 1) * BUFMIN;
    memory->grow(recvcells, recv_buff_size, "kinetics::recvcells");
  }
  if (is_intersection_valid(recvgrid)) {
    add_cells(basegrid, recvgrid, recvcells, nrecv);
    nrecv += n;
  }
  // identify which cells we need to send
  Grid sendgrid = intersect(extend(other), grid);
  n = cell_count(sendgrid);
  if (send_buff_size < nsend + n) { // not enough memory to fit cells
    send_buff_size += (n / BUFMIN + 1) * BUFMIN;
    memory->grow(sendcells, send_buff_size, "kinetics::sendcells");
  }
  if (is_intersection_valid(sendgrid)) {
    add_cells(basegrid, sendgrid, sendcells, nsend);
    nsend += n;
  }
}

/* ----------------------------------------------------------------------
   reset nutrient reaction array
 ------------------------------------------------------------------------- */

void FixKinetics::reset_nuR() {
  for (int k = 1; k < nnus+1; k++) {
    for (int j = 0; j < bgrids; j++) {
        nuR[k][j] = 0;
     }
   }
}

/* ----------------------------------------------------------------------
   reset convergence status
 ------------------------------------------------------------------------- */

void FixKinetics::reset_isConv() {
  for (int i = 1; i <= nnus; i++) {
    if (strcmp(bio->nuName[i], "h2o") == 0)
      nuConv[i] = true;
    else if (strcmp(bio->nuName[i], "h") == 0)
      nuConv[i] = true;
    else if (bio->diffCoeff[i] == 0)
      nuConv[i] = true;
    else
      nuConv[i] = false;
  }
}

/* ---------------------------------------------------------------------- */

int FixKinetics::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"demflag") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal fix_modify command");
    if (strcmp(arg[1],"no") == 0) demflag = 0;
    else if (strcmp(arg[1],"yes") == 0) demflag = 1;
    else error->all(FLERR,"Illegal demflag parameter (yes or no)");
    return 2;
  }
  return 0;
}

#ifdef OUTPUT_GRID
void FixKinetics::output_nutrient_info(int nutrient) {
  if (nutrient < 1 || nutrient > nnus)
    return;

  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(subn[0] + 1, subn[1] + 1, subn[2] + 1);
  image->SetOrigin(sublo[0], sublo[1], sublo[2]);
  image->SetSpacing(stepx, stepy, stepz);

  vtkSmartPointer<vtkDoubleArray> concentration = vtkSmartPointer<vtkDoubleArray>::New();
  concentration->SetName("Concentration");
  concentration->SetNumberOfComponents(1);

  vtkSmartPointer<vtkDoubleArray> uptake = vtkSmartPointer<vtkDoubleArray>::New();
  uptake->SetName("Uptake");
  uptake->SetNumberOfComponents(1);
  
  vtkSmartPointer<vtkDoubleArray> gas;
  vtkSmartPointer<vtkDoubleArray> act;

  if (energy != NULL) {
    gas = vtkSmartPointer<vtkDoubleArray>::New();
    gas->SetName("Gas");
    gas->SetNumberOfComponents(1);

    act = vtkSmartPointer<vtkDoubleArray>::New();
    act->SetName("Activity");
    act->SetNumberOfComponents(5);
  }

  int index = 0;
  for (int z = 0; z < subn[2]; z++) {
    for (int y = 0; y < subn[1]; y++) {
      for (int x = 0; x < subn[0]; x++, index++) {
	concentration->InsertNextTuple1(nuS[nutrient][index]);
	uptake->InsertNextTuple1(nuR[nutrient][index]);
	if (energy != NULL) {
	  gas->InsertNextTuple1(qGas[nutrient][index]);
	  double tuple[5];
	  for (int i = 0; i < 5; i++)
	    tuple[i] = activity[nutrient][i][index];
	  act->InsertNextTuple(tuple);
	}
      }
    }
  }

  image->GetCellData()->AddArray(concentration);
  image->GetCellData()->AddArray(uptake);
  if (energy != NULL) {
    image->GetCellData()->AddArray(gas);
    image->GetCellData()->AddArray(act);
  }

  const int maxsize = 256; 
  char *filename = new char[maxsize];
  sprintf(filename, "%s_%d_%ld.vti", bio->nuName[nutrient], comm->me, update->ntimestep);

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(filename);
  writer->SetInputData(image);
  writer->Write();

  delete [] filename;
}

void FixKinetics::output_bacteria_info(int type) {
  if (type < 1 || type > atom->ntypes)
    return;

  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(subn[0] + 1, subn[1] + 1, subn[2] + 1);
  image->SetOrigin(sublo[0], sublo[1], sublo[2]);
  image->SetSpacing(stepx, stepy, stepz);

  vtkSmartPointer<vtkDoubleArray> yield = vtkSmartPointer<vtkDoubleArray>::New();
  yield->SetName("Yield");
  yield->SetNumberOfComponents(1);

  vtkSmartPointer<vtkDoubleArray> cat = vtkSmartPointer<vtkDoubleArray>::New();
  cat->SetName("Uptake");
  cat->SetNumberOfComponents(1);

  vtkSmartPointer<vtkDoubleArray> ana = vtkSmartPointer<vtkDoubleArray>::New();
  ana->SetName("Uptake");
  ana->SetNumberOfComponents(1);

  int index = 0;
  for (int z = 0; z < subn[2]; z++) {
    for (int y = 0; y < subn[1]; y++) {
      for (int x = 0; x < subn[0]; x++, index++) {
	yield->InsertNextTuple1(gYield[type][index]);
	cat->InsertNextTuple1(DRGCat[type][index]);
	ana->InsertNextTuple1(DRGAn[type][index]);
      }
    }
  }

  image->GetCellData()->AddArray(yield);
  image->GetCellData()->AddArray(cat);
  image->GetCellData()->AddArray(ana);

  const int maxsize = 256; 
  char *filename = new char[maxsize];
  sprintf(filename, "%s_%d_%ld.vti", bio->typeName[type], comm->me, update->ntimestep);

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(filename);
  writer->SetInputData(image);
  writer->Write();

  delete [] filename;
}
#endif // OUTPUT_GRID
