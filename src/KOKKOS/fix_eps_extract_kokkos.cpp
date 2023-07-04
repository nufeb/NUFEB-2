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

#include "fix_eps_extract_kokkos.h"
#include "atom_kokkos.h"
#include "atom_vec_kokkos.h"
#include "comm.h"
#include "force.h"
#include "error.h"
#include "domain.h"
#include "group.h"
#include "lammps.h"
#include "compute.h"
#include "modify.h"
#include "kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixEPSExtractKokkos<DeviceType>::FixEPSExtractKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixEPSSecretion(lmp, narg, arg), rand_pool(seed + comm->me)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | RADIUS_MASK |
      OUTER_MASS_MASK | OUTER_RADIUS_MASK | OMEGA_MASK | TORQUE_MASK | BIOMASS_MASK | TAG_MASK;
  datamask_modify = datamask_read;
  divide_list = nullptr;
  nlocal = atomKK->nlocal;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixEPSExtractKokkos<DeviceType>::~FixEPSExtractKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_divide_list,divide_list);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixEPSExtractKokkos<DeviceType>::init()
{
  grow_arrays(atomKK->nlocal);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixEPSExtractKokkos<DeviceType>::grow_arrays(int nlocal)
{
  k_divide_list.template sync<LMPHostType>();
  memoryKK->grow_kokkos(k_divide_list, divide_list, nlocal, "nufeb/eps_extract:divide_list");
  d_divide_list = k_divide_list.template view<DeviceType>();
  h_divide_list = k_divide_list.template view<LMPHostType>();
  k_divide_list.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixEPSExtractKokkos<DeviceType>::compute()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->sync(Host, MASK_MASK | RADIUS_MASK | OUTER_RADIUS_MASK);

  for (int i = 0; i < 3; i++) {
    boxlo[i] = domain->boxlo[i];
    boxhi[i] = domain->boxhi[i];
  }

  if (nlocal < atomKK->nlocal) {
    grow_arrays(atomKK->nlocal);
    nlocal = atomKK->nlocal;
  }

  d_x = atomKK->k_x.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_tag = atomKK->k_tag.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_omega = atomKK->k_omega.view<DeviceType>();
  d_torque = atomKK->k_torque.view<DeviceType>();

  d_rmass = atomKK->k_rmass.template view<DeviceType>();
  d_biomass = atomKK->k_biomass.template view<DeviceType>();
  d_radius = atomKK->k_radius.template view<DeviceType>();
  d_outer_mass = atomKK->k_outer_mass.template view<DeviceType>();
  d_outer_radius = atomKK->k_outer_radius.template view<DeviceType>();

  h_mask = atomKK->k_mask.view<LMPHostType>();
  h_radius = atomKK->k_radius.view<LMPHostType>();
  h_outer_radius = atomKK->k_outer_radius.view<LMPHostType>();

  // create new atom on host
  for (int i = 0; i < nlocal; i++) {
    if (h_mask(i) & groupbit) {
      if (h_outer_radius(i) / h_radius(i) > ratio) {
	double *coord = new double[3];
	coord[0] = coord[1] = coord[2] = 0.0;
	atomKK->avec->create_atom(type, coord);
	h_divide_list(i) = atomKK->nlocal - 1;
	delete[] coord;
      } else {
	h_divide_list(i) = -1;
      }
    }
  }

  k_divide_list.template modify<LMPHostType>();
  k_divide_list.template sync<DeviceType>();
  eps_mask = group->bitmask[ieps];

  // update atom attributes on device
  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixEPSExtractComputeTag>(0, nlocal), f);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);
  atomKK->sync(Host,TAG_MASK);

  bigint nblocal = atomKK->nlocal;
  MPI_Allreduce(&nblocal, &atomKK->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atomKK->natoms < 0 || atomKK->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atomKK->tag_enable) atomKK->tag_extend();
  atomKK->tag_check();

  if (atomKK->map_style) {
    atomKK->nghost = 0;
    atomKK->map_init();
    atomKK->map_set();
  }

  atomKK->modified(Host,TAG_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixEPSExtractKokkos<DeviceType>::Functor::Functor(FixEPSExtractKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit), eps_density(ptr->eps_density), d_divide_list(ptr->d_divide_list),
  d_x(ptr->d_x), d_f(ptr->d_f), d_v(ptr->d_v), d_mask(ptr->d_mask), d_tag(ptr->d_tag),
  d_omega(ptr->d_omega), d_torque(ptr->d_torque), d_rmass(ptr->d_rmass), d_biomass(ptr->d_biomass),
  d_radius(ptr->d_radius), d_outer_mass(ptr->d_outer_mass), d_outer_radius(ptr->d_outer_radius),
  rand_pool(ptr->rand_pool), eps_mask(ptr->eps_mask)
{
  for (int i = 0; i < 3; i++) {
    boxlo[i] = ptr->boxlo[i];
    boxhi[i] = ptr->boxhi[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixEPSExtractKokkos<DeviceType>::Functor::operator()(FixEPSExtractComputeTag, int i) const
{
  if (d_mask(i) & groupbit) {
    if (d_divide_list(i) > 0) {
      rand_type rand_gen = rand_pool.get_state();
      static const double MY_PI  = 3.14159265358979323846; // pi
      double DELTA = 1.005;

      d_outer_mass(i) = (4.0 * MY_PI / 3.0) *
	  ((d_outer_radius(i) * d_outer_radius(i) * d_outer_radius(i)) -
	      (d_radius(i) * d_radius(i) * d_radius(i))) * eps_density;

      double split = 0.4 + (rand_gen.drand() * 0.2); //(random->uniform() * 0.2);

      double new_outer_mass = d_outer_mass(i) * split;
      double eps_mass = d_outer_mass(i) - new_outer_mass;

      d_outer_mass(i) = new_outer_mass;

      double density = d_rmass(i) / (4.0 * MY_PI / 3.0 * d_radius(i) * d_radius(i) * d_radius(i));
      d_outer_radius(i) = pow((3.0 / (4.0 * MY_PI)) * ((d_rmass(i) / density) + (d_outer_mass(i) / eps_density)), (1.0 / 3.0));

      double theta = rand_gen.drand() * 2 * MY_PI; //random->uniform() * 2 * MY_PI;
      double phi = rand_gen.drand() * MY_PI; // random->uniform() * (MY_PI);

      double oldx = d_x(i,0);
      double oldy = d_x(i,1);
      double oldz = d_x(i,2);

      // update attributes of the new EPS atom
      double child_radius = pow(((6 * eps_mass) / (eps_density * MY_PI)), (1.0 / 3.0)) * 0.5;
      double newx = oldx - (child_radius * cos(theta) * sin(phi) * DELTA);
      double newy = oldy - (child_radius * sin(theta) * sin(phi) * DELTA);
      double newz = oldz - (child_radius * cos(phi) * DELTA);

      if (newx - child_radius < boxlo[0]) {
        newx = boxlo[0] + child_radius;
      } else if (newx + d_outer_radius(i) > boxhi[0]) {
        newx = boxhi[0] - child_radius;
      }
      if (newy - child_radius < boxlo[1]) {
        newy = boxlo[1] + child_radius;
      } else if (newy + d_outer_radius(i) > boxhi[1]) {
        newy = boxhi[1] - child_radius;
      }
      if (newz - child_radius < boxlo[2]) {
        newz = boxlo[2] + child_radius;
      } else if (newz + d_outer_radius(i) > boxhi[2]) {
        newz = boxhi[2] - child_radius;
      }

      int j = d_divide_list[i];

      d_x(j,0) = newx;
      d_x(j,1) = newy;
      d_x(j,2) = newz;

      d_tag(j) = 0;
      d_mask(j) = 1 | eps_mask;
      d_v(j,0) = d_v(i,0);
      d_v(j,1) = d_v(i,1);
      d_v(j,2) = d_v(i,2);
      d_f(j,0) = d_f(i,0);
      d_f(j,1) = d_f(i,1);
      d_f(j,2) = d_f(i,2);
      d_omega(j,0) = d_omega(i,0);
      d_omega(j,1) = d_omega(i,1);
      d_omega(j,2) = d_omega(i,2);
      d_torque(j,0) = d_torque(i,0);
      d_torque(j,1) = d_torque(i,1);
      d_torque(j,2) = d_torque(i,2);
      d_rmass(j) = eps_mass;
      d_biomass(j) = 1.0;
      d_radius(j) = child_radius;
      d_outer_mass(j) = 0.0;
      d_outer_radius(j) = child_radius;

      rand_pool.free_state(rand_gen);
    }
  }
}


/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixEPSExtractKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixEPSExtractKokkos<LMPHostType>;
#endif
}
