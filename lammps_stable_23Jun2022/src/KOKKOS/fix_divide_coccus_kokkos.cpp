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

#include "fix_divide_coccus_kokkos.h"
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
FixDivideCoccusKokkos<DeviceType>::FixDivideCoccusKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixDivideCoccus(lmp, narg, arg),rand_pool(seed + comm->me)
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
FixDivideCoccusKokkos<DeviceType>::~FixDivideCoccusKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_divide_list,divide_list);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDivideCoccusKokkos<DeviceType>::init()
{
  grow_arrays(atomKK->nlocal);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDivideCoccusKokkos<DeviceType>::grow_arrays(int nlocal)
{
  k_divide_list.template sync<LMPHostType>();
  memoryKK->grow_kokkos(k_divide_list, divide_list, nlocal, "nufeb/division/coccus:divide_list");
  d_divide_list = k_divide_list.template view<DeviceType>();
  h_divide_list = k_divide_list.template view<LMPHostType>();
  k_divide_list.template modify<LMPHostType>();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixDivideCoccusKokkos<DeviceType>::compute()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->sync(Host, MASK_MASK | RADIUS_MASK | TYPE_MASK);

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
  h_type = atomKK->k_type.view<LMPHostType>();
  h_radius = atomKK->k_radius.view<LMPHostType>();

  // create new atom on host
  for (int i = 0; i < nlocal; i++) {
    if (h_mask(i) & groupbit) {
      if (h_radius(i) * 2 >= diameter) {
	double *coord = new double[3];
	coord[0] = coord[1] = coord[2] = 0.0;
	atomKK->avec->create_atom(h_type(i), coord);
	h_divide_list(i) = atomKK->nlocal - 1;
	delete[] coord;
      } else {
	h_divide_list(i) = -1;
      }
    }
  }

  k_divide_list.template modify<LMPHostType>();
  k_divide_list.template sync<DeviceType>();

  // update atom attributes on device
  copymode = 1;
  Functor f(this);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<
    DeviceType,
    FixDivideCoccusComputeTag>(0, nlocal), f);
  copymode = 0;

  atomKK->modified(execution_space,datamask_modify);

  for (int i = 0; i < nlocal; i++) {
    int j = h_divide_list(i);
    if (j > 0) {
      modify->create_attribute(j);

      for (int m = 0; m < modify->nfix; m++)
	modify->fix[m]->update_arrays(i, j);
    }
  }

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
FixDivideCoccusKokkos<DeviceType>::Functor::Functor(FixDivideCoccusKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit), eps_density(ptr->eps_density), d_divide_list(ptr->d_divide_list),
  d_x(ptr->d_x), d_f(ptr->d_f), d_v(ptr->d_v), d_mask(ptr->d_mask), d_tag(ptr->d_tag),
  d_omega(ptr->d_omega), d_torque(ptr->d_torque), d_rmass(ptr->d_rmass), d_biomass(ptr->d_biomass),
  d_radius(ptr->d_radius), d_outer_mass(ptr->d_outer_mass), d_outer_radius(ptr->d_outer_radius),
  rand_pool(ptr->rand_pool)
{
  for (int i = 0; i < 3; i++) {
    boxlo[i] = ptr->boxlo[i];
    boxhi[i] = ptr->boxhi[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixDivideCoccusKokkos<DeviceType>::Functor::operator()(FixDivideCoccusComputeTag, int i) const
{
  if (d_mask(i) & groupbit) {
    if (d_divide_list(i) > 0) {
      rand_type rand_gen = rand_pool.get_state();
      static const double MY_PI  = 3.14159265358979323846; // pi
      double DELTA = 1.005;

      double density = d_rmass(i) /
	  (4.0 * MY_PI / 3.0 * d_radius(i) * d_radius(i) * d_radius(i));
      double split = 0.4 + (rand_gen.drand() * 0.2); //(random->uniform() * 0.2);
      double imass = d_rmass(i) * split;
      double jmass = d_rmass(i) - imass;

      double iouter_mass = d_outer_mass(i) * split;
      double jouter_mass = d_outer_mass(i) - iouter_mass;

      double theta = rand_gen.drand() * 2 * MY_PI; //random->uniform() * 2 * MY_PI;
      double phi = rand_gen.drand() * MY_PI; // random->uniform() * (MY_PI);

      double oldx = d_x(i,0);
      double oldy = d_x(i,1);
      double oldz = d_x(i,2);

      // update daughter cell i
      d_rmass(i) = imass;
      d_outer_mass(i) = iouter_mass;
      d_radius(i) = pow(((6 * d_rmass(i)) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
      d_outer_radius(i) = pow((3.0 / (4.0 * MY_PI)) * ((d_rmass(i) / density) + (iouter_mass / eps_density)), (1.0 / 3.0));

      double newx = oldx + (d_radius(i) * cos(theta) * sin(phi) * DELTA);
      double newy = oldy + (d_radius(i) * sin(theta) * sin(phi) * DELTA);
      double newz = oldz + (d_radius(i) * cos(phi) * DELTA);

      if (newx - d_radius(i) < boxlo[0]) {
        newx = boxlo[0] + d_radius(i);
      } else if (newx + d_radius(i) > boxhi[0]) {
        newx = boxhi[0] - d_radius(i);
      }
      if (newy - d_radius(i) < boxlo[1]) {
        newy = boxlo[1] + d_radius(i);
      } else if (newy + d_radius(i) > boxhi[1]) {
        newy = boxhi[1] - d_radius(i);
      }
      if (newz - d_radius(i) < boxlo[2]) {
        newz = boxlo[2] + d_radius(i);
      } else if (newz + d_radius(i) > boxhi[2]) {
        newz = boxhi[2] - d_radius(i);
      }

      d_x(i,0) = newx;
      d_x(i,1) = newy;
      d_x(i,2) = newz;

      int j = d_divide_list[i];

      // create daughter cell j
      double jradius = pow(((6 * jmass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
      double jouter_radius = pow((3.0 / (4.0 * MY_PI)) * ((jmass / density) + (jouter_mass / eps_density)), (1.0 / 3.0));
      newx = oldx - (jradius * cos(theta) * sin(phi) * DELTA);
      newy = oldy - (jradius * sin(theta) * sin(phi) * DELTA);
      newz = oldz - (jradius * cos(phi) * DELTA);

      if (newx - jradius < boxlo[0]) {
        newx = boxlo[0] + jradius;
      } else if (newx + jradius > boxhi[0]) {
        newx = boxhi[0] - jradius;
      }
      if (newy - jradius < boxlo[1]) {
        newy = boxlo[1] + jradius;
      } else if (newy + jradius > boxhi[1]) {
        newy = boxhi[1] - jradius;
      }
      if (newz - jradius < boxlo[2]) {
        newz = boxlo[2] + jradius;
      } else if (newz + jradius > boxhi[2]) {
        newz = boxhi[2] - jradius;
      }

      d_x(j,0) = newx;
      d_x(j,1) = newy;
      d_x(j,2) = newz;

      d_tag(j) = 0;
      d_mask(j) = d_mask(i);
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

      d_rmass(j) = jmass;
      d_biomass(j) = d_biomass(i);
      d_radius(j) = jradius;
      d_outer_mass(j) = jouter_mass;
      d_outer_radius(j) = jouter_radius;

      rand_pool.free_state(rand_gen);
    }
  }
}


/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixDivideCoccusKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixDivideCoccusKokkos<LMPHostType>;
#endif
}
