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

#include <math.h>
#include <string.h>
#include "fix_eps_adhesion_kokkos.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "neigh_list_kokkos.h"
#include "pair_kokkos.h"
#include "force.h"
#include "pair.h"
#include "group.h"
#include "lammps.h"
#include "kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{DEFAULT,SQUARE};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixEPSAdhesionKokkos<DeviceType>::FixEPSAdhesionKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEPSAdhesion(lmp, narg, arg)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  virial_flag = 0;

  datamask_read = X_MASK | F_MASK | MASK_MASK | RMASS_MASK | RADIUS_MASK | OUTER_MASS_MASK | OUTER_RADIUS_MASK;
  datamask_modify = F_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixEPSAdhesionKokkos<DeviceType>::post_force(int)
{
  int inum = force->pair->list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(force->pair->list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  d_x = atomKK->k_x.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();
  d_outer_mass = atomKK->k_outer_mass.view<DeviceType>();
  d_outer_radius = atomKK->k_outer_radius.view<DeviceType>();

  epsmask = group->bitmask[ieps];
  nlocal = atom->nlocal;

  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  copymode = 1;
  Functor f(this);
  if (lmp->kokkos->neighflag == HALF) {
    if (force->newton_pair) {
      if (disp == DEFAULT) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALF, 0, 0> >(0, inum), f);
      } else if (disp == SQUARE) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALF, 0, 1> >(0, inum), f);
      }
    } else {
      if (disp == DEFAULT) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALF, 1, 0> >(0, inum), f);
      } else if (disp == SQUARE) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALF, 1, 1> >(0, inum), f);
      }
    }
  } else if (lmp->kokkos->neighflag == HALFTHREAD) {
    if (force->newton_pair) {
      if (disp == DEFAULT) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALFTHREAD, 0, 0> >(0, inum), f);
      } else if (disp == SQUARE) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALFTHREAD, 0, 1> >(0, inum), f);
      }
    } else {
      if (disp == DEFAULT) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALFTHREAD, 1, 0> >(0, inum), f);
      } else if (disp == SQUARE) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<HALFTHREAD, 1, 1> >(0, inum), f);
      }
    }
  } else { // lmp->kokkos->neighflag == FULL
    if (force->newton_pair) {
      if (disp == DEFAULT) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<FULL, 0, 0> >(0, inum), f);
      } else if (disp == SQUARE) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<FULL, 0, 1> >(0, inum), f);
      }
    } else {
      if (disp == DEFAULT) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<FULL, 1, 0> >(0, inum), f);
      } else if (disp == SQUARE) {
	Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixEPSAdhesionTag<FULL, 1, 1> >(0, inum), f);
      }
    }
  }
  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixEPSAdhesionKokkos<DeviceType>::Functor::Functor(FixEPSAdhesionKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit), epsmask(ptr->epsmask), nlocal(ptr->nlocal), ke(ptr->ke),
  d_neighbors(ptr->d_neighbors), d_ilist(ptr->d_ilist), d_numneigh(ptr->d_numneigh),
  d_x(ptr->d_x), d_f(ptr->d_f), d_mask(ptr->d_mask), d_rmass(ptr->d_rmass),
  d_radius(ptr->d_radius), d_outer_mass(ptr->d_outer_mass),
  d_outer_radius(ptr->d_outer_radius) {}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int NEIGHFLAG, int NEWTON_PAIR, int DISP>
KOKKOS_INLINE_FUNCTION
void FixEPSAdhesionKokkos<DeviceType>::Functor::operator()(FixEPSAdhesionTag<NEIGHFLAG, NEWTON_PAIR, DISP>, int ii) const
{
  // The f array is atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = d_f;
  
  const int i = d_ilist[ii];
  const int jnum = d_numneigh[i];

  if (d_mask[i] & groupbit) {
    double xtmp = d_x(i,0);
    double ytmp = d_x(i,1);
    double ztmp = d_x(i,2);

    double radi = d_radius[i];
    double oradi = d_outer_radius[i];

    double epsi = 0;
    if (d_mask[i] & epsmask)
      epsi = d_rmass[i];
    else
      epsi = d_outer_mass[i];

    double fx = 0.0; 
    double fy = 0.0; 
    double fz = 0.0; 
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      double delx = xtmp - d_x(j,0);
      double dely = ytmp - d_x(j,1);
      double delz = ztmp - d_x(j,2);
      double rsq = delx * delx + dely * dely + delz * delz;

      double radj = d_radius[j];
      double oradj = d_outer_radius[j];

      double epsj = 0;
      if (d_mask[j] & epsmask)
	epsj = d_rmass[j];
      else
	epsj = d_outer_mass[j];

      double radsum = radi + radj;
      double oradsum = oradi + oradj;
      double masssum = epsi + epsj;
      double r = sqrt(rsq);
      // double del = r - radsum;
      double del = r - 0.5 * (radsum + oradsum);
      double rinv = 1 / r;

      double ccel = 0;
      if (r > radsum && r < oradsum) {
	if (DISP == DEFAULT)
	  ccel = -masssum * ke * del;
	if (DISP == SQUARE)
	  ccel = -masssum * ke * (radsum / r) * (radsum / r);
      }

      double ccelx = delx * ccel * rinv;
      double ccely = dely * ccel * rinv;
      double ccelz = delz * ccel * rinv;

      fx += ccelx;
      fy += ccely;
      fz += ccelz;

      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
	a_f(j,0) -= ccelx;
	a_f(j,1) -= ccely;
	a_f(j,2) -= ccelz;
      }
    }
    
    a_f(i,0) += fx;
    a_f(i,1) += fy;
    a_f(i,2) += fz;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixEPSAdhesionKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixEPSAdhesionKokkos<LMPHostType>;
#endif
}
