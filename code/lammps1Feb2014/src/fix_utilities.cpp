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
#include "stdlib.h"
#include "fix_utilities.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "math.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "mpi.h"
#include "comm.h"
#include "memory.h"
#include "input.h"
#include "variable.h"
#include <iostream>
#include "math_const.h"
#include <vector>


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using namespace std;

#define DELTA 10000;
/* ---------------------------------------------------------------------- */

FixUtilities::FixUtilities(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal fix utilities command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix utilities command: calling steps should be positive integer");
}

FixUtilities::~FixUtilities()
{
}

/* ---------------------------------------------------------------------- */

int FixUtilities::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixUtilities::init()
{
//	int irequest = neighbor->request((void *) this);
//  neighbor->requests[irequest]->pair = 0;
//  neighbor->requests[irequest]->fix = 1;

}

///* ---------------------------------------------------------------------- */
//
//void FixUtilities::init_list(int id, NeighList *ptr)
//{
//  list = ptr;
//}


void FixUtilities::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  list.clear();
  int nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;
  fourThirdsPI = 4.0*MY_PI/3.0;

  visit = new int[nall]();
  vector<double> floc_vol;
  vector<int> floc_sizes;

  neighbor_list();
  
  for (int i = 0; i < nall; i++) {
    if (visit[i] == 0) {
      int size = 1;
      double radius = atom->radius[i];
      double volume = fourThirdsPI * radius * radius * radius;
      double vol = volume;

      floc_size (vol, i, size);
      floc_vol.push_back(vol);
      floc_sizes.push_back(size);
    }
  }

  int n = 0;
  for (auto i = floc_vol.begin(); i != floc_vol.end(); ++i) {
      n++;
      cout << "Floc " << n << " size = " << floc_sizes.at(n-1) << endl;
      cout <<  *i << ' ';
      cout << endl;
  }

  delete[] visit;
}

void FixUtilities::floc_size (double &vol, int bac, int &size) {
  visit[bac] = 1;

  for (int const& j: list.at(bac)) {
//    int j = jlist[jj];
    if (visit[j] == 0) {
      double radius = atom->radius[j];
      double volume = fourThirdsPI * radius * radius * radius;
      vol = vol + volume;
      size ++;
      floc_size (vol, j, size);
    }
  }
}

void FixUtilities::neighbor_list () {

  for(int i = 0; i < nall; i++){
    vector<int> subList;
    for(int j = 0; j < nall; j++){
      if(i != j) {
        double xd = atom->x[i][0] - atom->x[j][0];
        double yd = atom->x[i][1] - atom->x[j][1];
        double zd = atom->x[i][2] - atom->x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (atom->radius[i] + atom->radius[j] + 1.0e-6) * (atom->radius[i] + atom->radius[j]+ 1.0e-6);

        if (rsq <= cut) subList.push_back(j);
      }
    }
    list.push_back(subList);
  }
}
