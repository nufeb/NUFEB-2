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

#include <cstdio>
#include <cstring>
#include "fix_nve_limit_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVELimitKokkos<DeviceType>::FixNVELimitKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVELimit(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = X_MASK | V_MASK;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNVELimitKokkos<DeviceType>::initial_integrate(int)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  Functor f(this);
  if (atom->rmass) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixNVELimitInitialRMassTag>(0, nlocal), f);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixNVELimitInitialTag>(0, nlocal), f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVELimitKokkos<DeviceType>::Functor::Functor(FixNVELimitKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit), dtf(ptr->dtf), dtv(ptr->dtv), vlimitsq(ptr->vlimitsq),
  x(ptr->x), v(ptr->v), f(ptr->f), rmass(ptr->rmass), mass(ptr->mass),
  type(ptr->type), mask(ptr->mask) {}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVELimitKokkos<DeviceType>::Functor::operator()(FixNVELimitInitialTag, int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / mass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    double vsq = v(i, 0)*v(i, 0) + v(i, 1)*v(i, 1) + v(i, 2)*v(i, 2);
    if (vsq > vlimitsq) {
      // ncount++;
      double scale = sqrt(vlimitsq/vsq);
      v(i, 0) *= scale;
      v(i, 1) *= scale;
      v(i, 2) *= scale;
    }
    
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVELimitKokkos<DeviceType>::Functor::operator()(FixNVELimitInitialRMassTag, int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / rmass[i];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    double vsq = v(i, 0)*v(i, 0) + v(i, 1)*v(i, 1) + v(i, 2)*v(i, 2);
    if (vsq > vlimitsq) {
      // ncount++;
      double scale = sqrt(vlimitsq/vsq);
      v(i, 0) *= scale;
      v(i, 1) *= scale;
      v(i, 2) *= scale;
    }

    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVELimitKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space,V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  atomKK->modified(execution_space,V_MASK);

  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  Functor f(this);
  if (atom->rmass) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixNVELimitFinalRMassTag>(0, nlocal), f);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixNVELimitFinalTag>(0, nlocal), f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVELimitKokkos<DeviceType>::Functor::operator()(FixNVELimitFinalTag, int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / mass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    double vsq = v(i, 0)*v(i, 0) + v(i, 1)*v(i, 1) + v(i, 2)*v(i, 2);
    if (vsq > vlimitsq) {
      // ncount++;
      double scale = sqrt(vlimitsq/vsq);
      v(i, 0) *= scale;
      v(i, 1) *= scale;
      v(i, 2) *= scale;
    }
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVELimitKokkos<DeviceType>::Functor::operator()(FixNVELimitFinalRMassTag, int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / rmass[i];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    double vsq = v(i, 0)*v(i, 0) + v(i, 1)*v(i, 1) + v(i, 2)*v(i, 2);
    if (vsq > vlimitsq) {
      // ncount++;
      double scale = sqrt(vlimitsq/vsq);
      v(i, 0) *= scale;
      v(i, 1) *= scale;
      v(i, 2) *= scale;
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double FixNVELimitKokkos<DeviceType>::compute_scalar()
{
  error->all(FLERR, "fix nve/limit/kk has not yet implemented compute_scalar.");
  return 0.0;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixNVELimitKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixNVELimitKokkos<LMPHostType>;
#endif
}

