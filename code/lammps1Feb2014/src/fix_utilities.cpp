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
}

void FixUtilities::post_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;

  list.clear();
  int nlocal = atom->nlocal;
  nall = nlocal + atom->nghost;
  fourThirdsPI = 4.0*MY_PI/3.0;

  visit = new int[nall]();
  vector <vector<int>> flocs;
  int *id = new int[nall]();
  vector<int> parent_floc;
  int max_size = 0;

  //build neighbor list
  neighbor_list();
  
  for (int i = 0; i < nall; i++) {
    if (visit[i] == 0) {
      vector<int> floc;
      floc.push_back(i);
      get_floc (i, floc);

      //get parent floc
      if (floc.size() > max_size) {
        max_size = floc.size();
        parent_floc = floc;
      }
    }
  }

  //remove parent floc
  flocs.erase(remove(flocs.begin(), flocs.end(), parent_floc), flocs.end());

  for (const vector<int> &floc : flocs) {
    //new detached floc
    if (id[floc.at(0)] == 0) {
      double vol;
      double x, y ,z;
      printf("New floc at nstep = %i:\n", update->nsteps);
      for (const int &i : floc) {
        if (id[0] == 1) cout << "warning: floc is mixed with pre-detached atoms" << endl;
        id[i] = 1;
        double radius = atom->radius[i];
        double volume = fourThirdsPI * radius * radius * radius;
        vol += volume;
        x += atom->x[i][0];
        y += atom->x[i][1];
        z += atom->x[i][2];
      }
      double size = floc.size();
      x = x/size;
      y = y/size;
      z = z/size;
      printf("  Size = %e:\n", size);
      printf("  Average position x = %e, y = %e, z = %e:\n", x, y, z);
      printf("  Volume = %e:\n", vol);
    }
  }

  delete[] visit;
  delete[] id;
}

void FixUtilities::get_floc (int bac, vector<int>& floc) {
  visit[bac] = 1;

  for (int const& j: list.at(bac)) {

    if (visit[j] == 0) {
      floc.push_back(j);
      get_floc (j, floc);
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
