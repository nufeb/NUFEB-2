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

#include "fix_wall_gran_kokkos.h"
#include "atom_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "atom_vec_kokkos.h"
#include "atom_masks.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER,REGION};
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY,BONDED_HISTORY};
enum{NONE,CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallGranKokkos<DeviceType>::FixWallGranKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | OMEGA_MASK | TORQUE_MASK | RADIUS_MASK | RMASS_MASK | MASK_MASK;
  datamask_modify = F_MASK | TORQUE_MASK;

  memory->destroy(history_one);
  history_one = NULL;
  grow_arrays(atom->nmax);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallGranKokkos<DeviceType>::~FixWallGranKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_history_one, history_one);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::init()
{
  FixWallGran::init();

  if (fix_rigid)
    error->all(FLERR, "wall/gran/kk not yet compatible with rigid.");  
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::post_force(int /*vflag*/)
{
  // do not update shear history during setup

  history_update = 1;
  if (update->setupflag) history_update = 0;

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  wlo = lo;
  whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = lo + amplitude - amplitude*cos(arg);
      whi = hi + amplitude - amplitude*cos(arg);
    }
    vwall[axis] = amplitude*omega*sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  copymode = 1;

  d_x = atomKK->k_x.view<DeviceType>();
  d_v = atomKK->k_v.view<DeviceType>();
  d_omega = atomKK->k_omega.view<DeviceType>();
  d_f = atomKK->k_f.view<DeviceType>();
  d_torque = atomKK->k_torque.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  d_rmass = atomKK->k_rmass.view<DeviceType>();
  d_radius = atomKK->k_radius.view<DeviceType>();
  int nlocal = atom->nlocal;

  if (pairstyle == HOOKE)
    error->all(FLERR, "wall/gran/kk doesn't yet support hooke style.");
  else if (pairstyle == HOOKE_HISTORY) {
    Functor f(this);
    if (wallstyle == XPLANE) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixWallGranTag<XPLANE> >(0, nlocal), f);
    } else if (wallstyle == YPLANE) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixWallGranTag<YPLANE> >(0, nlocal), f);
    } else if (wallstyle == ZPLANE) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixWallGranTag<ZPLANE> >(0, nlocal), f);
    } else if (wallstyle == ZCYLINDER) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, FixWallGranTag<ZCYLINDER> >(0, nlocal), f);
    }
  }
  else if (pairstyle == HERTZ_HISTORY)
    error->all(FLERR, "wall/gran/kk doesn't yet support hertz/history style.");

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixWallGranKokkos<DeviceType>::Functor::Functor(FixWallGranKokkos<DeviceType> *ptr):
  groupbit(ptr->groupbit),
  history_update(ptr->history_update),
  use_history(ptr->use_history),
  axis(ptr->axis),
  wlo(ptr->wlo), whi(ptr->whi),
  kn(ptr->kn), kt(ptr->kt),
  gamman(ptr->gamman), gammat(ptr->gammat),
  xmu(ptr->xmu), dt(ptr->dt),
  vshear(ptr->vshear), wshear(ptr->wshear),
  cylradius(ptr->cylradius),
  d_x(ptr->d_x), d_v(ptr->d_v),
  d_omega(ptr->d_omega), d_f(ptr->d_f),
  d_torque(ptr->d_torque), d_mask(ptr->d_mask),
  d_rmass(ptr->d_rmass), d_radius(ptr->d_radius),
  d_history_one(ptr->d_history_one)
{
  for (int i = 0; i < 3; i++) {
    vwall[i] = ptr->vwall[i];
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <int WallStyle>
void FixWallGranKokkos<DeviceType>::Functor::operator()(FixWallGranTag<WallStyle>, int i) const
{
  double vwall_[3];
  vwall_[0] = vwall[0];
  vwall_[1] = vwall[1];
  vwall_[2] = vwall[2];
  
  if (d_mask[i] & groupbit) {
    X_FLOAT radius = d_radius(i);

    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    
    if (WallStyle == XPLANE) {
      X_FLOAT del1 = d_x(i,0) - wlo;
      double del2 = whi - d_x(i,0);
      if (del1 < del2) dx = del1;
      else dx = -del2;
    } else if (WallStyle == YPLANE) {
      double del1 = d_x(i,1) - wlo;
      double del2 = whi - d_x(i,1);
      if (del1 < del2) dy = del1;
      else dy = -del2;
    } else if (WallStyle == ZPLANE) {
      double del1 = d_x(i,2) - wlo;
      double del2 = whi - d_x(i,2);
      if (del1 < del2) dz = del1;
      else dz = -del2;
    } else if (WallStyle == ZCYLINDER) {
      double delxy = sqrt(d_x(i,0)*d_x(i,0) + d_x(i,1)*d_x(i,1));
      double delr = cylradius - delxy;
      if (delr > radius) {
    	dz = cylradius;
      } else {
    	dx = -delr/delxy * d_x(i,0);
    	dy = -delr/delxy * d_x(i,1);
     	if (wshear && axis != 2) {
    	  vwall_[0] += vshear * d_x(i,1)/delxy;
    	  vwall_[1] += -vshear * d_x(i,0)/delxy;
    	  vwall_[2] = 0.0;
    	}
      }
    }

    double rsq = dx*dx + dy*dy + dz*dz;

    if (rsq > radius*radius) {
      if (use_history)
    	for (int j = 0; j < 3; j++)
    	  d_history_one(i,j) = 0.0;
    } else {
      // meff = effective mass of sphere
      double meff = d_rmass(i);
      double r = sqrt(rsq);
      double rinv = 1.0/r;
      double rsqinv = 1.0/rsq;

      // relative translational velocity

      double vr1 = d_v(i,0) - vwall_[0];
      double vr2 = d_v(i,1) - vwall_[1];
      double vr3 = d_v(i,2) - vwall_[2];

      // normal component

      double vnnr = vr1*dx + vr2*dy + vr3*dz;
      double vn1 = dx*vnnr * rsqinv;
      double vn2 = dy*vnnr * rsqinv;
      double vn3 = dz*vnnr * rsqinv;

      // tangential component

      double vt1 = vr1 - vn1;
      double vt2 = vr2 - vn2;
      double vt3 = vr3 - vn3;

      // relative rotational velocity

      double wr1 = radius*d_omega(i,0) * rinv;
      double wr2 = radius*d_omega(i,1) * rinv;
      double wr3 = radius*d_omega(i,2) * rinv;

      // normal forces = Hookian contact + normal velocity damping

      double damp = meff*gamman*vnnr*rsqinv;
      double ccel = kn*(radius-r)*rinv - damp;

      // relative velocities

      double vtr1 = vt1 - (dz*wr2-dy*wr3);
      double vtr2 = vt2 - (dx*wr3-dz*wr1);
      double vtr3 = vt3 - (dy*wr1-dx*wr2);
      double vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
      vrel = sqrt(vrel);

      // shear history effects

      if (history_update) {
    	d_history_one(i,0) += vtr1*dt;
    	d_history_one(i,1) += vtr2*dt;
    	d_history_one(i,2) += vtr3*dt;
      }
      double shrmag = sqrt(d_history_one(i,0)*d_history_one(i,0) + d_history_one(i,1)*d_history_one(i,1) + d_history_one(i,2)*d_history_one(i,2));

      // rotate shear displacements

      double rsht = d_history_one(i,0)*dx + d_history_one(i,1)*dy + d_history_one(i,2)*dz;
      rsht = rsht*rsqinv;
      if (history_update) {
    	d_history_one(i,0) -= rsht*dx;
    	d_history_one(i,1) -= rsht*dy;
    	d_history_one(i,2) -= rsht*dz;
      }

      // tangential forces = shear + tangential velocity damping

      double fs1 = - (kt*d_history_one(i,0) + meff*gammat*vtr1);
      double fs2 = - (kt*d_history_one(i,1) + meff*gammat*vtr2);
      double fs3 = - (kt*d_history_one(i,2) + meff*gammat*vtr3);

      // rescale frictional displacements and forces if needed

      double fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
      double fn = xmu * fabs(ccel*r);

      if (fs > fn) {
    	if (shrmag != 0.0) {
    	  d_history_one(i,0) = (fn/fs) * (d_history_one(i,0) + meff*gammat*vtr1/kt) -
    	    meff*gammat*vtr1/kt;
    	  d_history_one(i,1) = (fn/fs) * (d_history_one(i,1) + meff*gammat*vtr2/kt) -
    	    meff*gammat*vtr2/kt;
    	  d_history_one(i,2) = (fn/fs) * (d_history_one(i,2) + meff*gammat*vtr3/kt) -
    	    meff*gammat*vtr3/kt;
    	  fs1 *= fn/fs ;
    	  fs2 *= fn/fs;
    	  fs3 *= fn/fs;
    	} else fs1 = fs2 = fs3 = 0.0;
      }

      // forces & torques

      double fx = dx*ccel + fs1;
      double fy = dy*ccel + fs2;
      double fz = dz*ccel + fs3;
      d_f(i,0) += fx;
      d_f(i,1) += fy;
      d_f(i,2) += fz;

      double tor1 = rinv * (dy*fs3 - dz*fs2);
      double tor2 = rinv * (dz*fs1 - dx*fs3);
      double tor3 = rinv * (dx*fs2 - dy*fs1);
      d_torque(i,0) -= radius*tor1;
      d_torque(i,1) -= radius*tor2;
      d_torque(i,2) -= radius*tor3;
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::grow_arrays(int nmax)
{
  if (use_history) {
    k_history_one.template sync<LMPHostType>(); // force reallocation on host 
    memoryKK->grow_kokkos(k_history_one,history_one,nmax,size_history,"wall/gran/kk:history_one");
    d_history_one = k_history_one.template view<DeviceType>();
    k_history_one.template modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::copy_arrays(int i, int j, int /*delflag*/)
{
  if (use_history) {
    k_history_one.template sync<LMPHostType>();
    for (int m = 0; m < size_history; m++)
      history_one[j][m] = history_one[i][m];
    k_history_one.template modify<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_exchange(int i, double *buf)
{
  k_history_one.template sync<LMPHostType>();

  int n = 0;
  for (int j = 0; j < size_history; j++)
    buf[n++] = history_one[i][j];
  return n;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int FixWallGranKokkos<DeviceType>::unpack_exchange(int nlocal, double *buf)
{
  int n = 0;
  for (int j = 0; j < size_history; j++)
    history_one[nlocal][j] = buf[n++];

  k_history_one.template modify<LMPHostType>();

  return n;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixWallGranKokkos_PackExchangeFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d_const _sendlist;
  typename AT::t_int_1d_const _copylist;
  typename AT::t_float_2d _history_one;
  typename AT::t_xfloat_1d_um _buf;
  const int _dnum;

  FixWallGranKokkos_PackExchangeFunctor(
    const typename AT::tdual_xfloat_2d &buf,
    const typename AT::tdual_int_1d &sendlist,
    const typename AT::tdual_int_1d &copylist,
    const typename AT::tdual_float_2d &history_one,
    const int &dnum):
    _sendlist(sendlist.template view<DeviceType>()),
    _copylist(copylist.template view<DeviceType>()),
    _history_one(history_one.template view<DeviceType>()),
    _dnum(dnum)
  {
    _buf = typename AT::t_xfloat_1d_um(buf.template view<DeviceType>().data(),buf.extent(0)*buf.extent(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &mysend) const {
    const int i = _sendlist(mysend);
    int m = mysend*_dnum;
    for (int v = 0; v < _dnum; v++) {
      _buf(m++) = _history_one(i,v);
    }
    const int j = _copylist(mysend);
    if (j > -1) {
      for (int v = 0; v < _dnum; v++) {
	_history_one(i,v) = _history_one(j,v);
      }
    }
  }
 };

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int FixWallGranKokkos<DeviceType>::pack_exchange_kokkos(
  const int &nsend,
  DAT::tdual_xfloat_2d &buf,
  DAT::tdual_int_1d k_sendlist,
  DAT::tdual_int_1d k_copylist,
  ExecutionSpace space, int dim,
  X_FLOAT lo, X_FLOAT hi)
{
  k_history_one.template sync<DeviceType>();
  Kokkos::parallel_for(
    nsend,
    FixWallGranKokkos_PackExchangeFunctor<DeviceType>(
      buf,k_sendlist,k_copylist,k_history_one,size_history));
  return nsend*size_history;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
struct FixWallGranKokkos_UnpackExchangeFunctor
{
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_xfloat_1d_um _buf;
  typename AT::t_float_2d _history_one;
  typename AT::t_int_1d _indices;
  const int _dnum;

  FixWallGranKokkos_UnpackExchangeFunctor(
    const typename AT::tdual_xfloat_2d buf,
    const typename AT::tdual_float_2d &history_one,
    const typename AT::tdual_int_1d &indices,
    const int &dnum):
    _history_one(history_one.template view<DeviceType>()),
    _indices(indices.template view<DeviceType>()),
    _dnum(dnum)
  {
    _buf = typename AT::t_xfloat_1d_um(buf.template view<DeviceType>().data(),buf.extent(0)*buf.extent(1));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    int index = _indices(i);
    if (index > 0) {
      int m = i*_dnum;
      for (int v = 0; v < _dnum; v++) {
	_history_one(index,v) = _buf(m++);
      }
    }
  }
};

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixWallGranKokkos<DeviceType>::unpack_exchange_kokkos(
  DAT::tdual_xfloat_2d &k_buf,
  DAT::tdual_int_1d &indices,int nrecv,
  int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
  ExecutionSpace space)
{
  Kokkos::parallel_for(
    nrecv/(atom->avec->size_border + atom->avec->size_velocity + 2),
    FixWallGranKokkos_UnpackExchangeFunctor<DeviceType>(
      k_buf,k_history_one,indices,size_history));

  k_history_one.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixWallGranKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class FixWallGranKokkos<LMPHostType>;
#endif
}
