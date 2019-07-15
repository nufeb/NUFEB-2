/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/kk,FixWallGranKokkos<LMPDeviceType>)
FixStyle(wall/gran/kk/device,FixWallGranKokkos<LMPDeviceType>)
FixStyle(wall/gran/kk/host,FixWallGranKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_WALL_GRAN_KOKKOS_H
#define LMP_FIX_WALL_GRAN_KOKKOS_H

#include "fix_wall_gran.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

template <int WallStyle>
struct FixWallGranTag {};
  
template<class DeviceType>
class FixWallGranKokkos : public FixWallGran, public KokkosBase {
 public:
  FixWallGranKokkos(class LAMMPS *, int, char **);
  ~FixWallGranKokkos();
  void init();
  void post_force(int);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_exchange_kokkos(const int &nsend,DAT::tdual_xfloat_2d &buf,
			   DAT::tdual_int_1d k_sendlist,
			   DAT::tdual_int_1d k_copylist,
			   ExecutionSpace space, int dim,
			   X_FLOAT lo, X_FLOAT hi);
  void unpack_exchange_kokkos(DAT::tdual_xfloat_2d &k_buf,
                              DAT::tdual_int_1d &indices,int nrecv,
                              int nlocal,int dim,X_FLOAT lo,X_FLOAT hi,
                              ExecutionSpace space);

  struct Functor
  {
    int groupbit;
    int history_update;
    int use_history;
    int axis;
    
    X_FLOAT wlo;
    X_FLOAT whi;

    double kn;
    double kt;
    double gamman;
    double gammat;
    double xmu;
    double dt;
    double vshear;
    double wshear;
    double cylradius;
    
    typedef ArrayTypes<DeviceType> AT;
    typename AT::t_x_array d_x;
    typename AT::t_v_array d_v;
    typename AT::t_v_array d_omega;
    typename AT::t_f_array d_f;
    typename AT::t_f_array d_torque;
    typename AT::t_int_1d d_mask;
    typename AT::t_float_1d d_rmass;
    typename AT::t_float_1d d_radius;
    typename AT::t_float_2d d_history_one;

    V_FLOAT vwall[3];
    
    Functor(FixWallGranKokkos *ptr);

    template <int WallStyle>
    KOKKOS_INLINE_FUNCTION
    void operator()(FixWallGranTag<WallStyle>, int) const;
  };
  
 protected:
  X_FLOAT wlo;
  X_FLOAT whi;
  V_FLOAT vwall[3];

  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_x_array d_x;
  typename AT::t_v_array d_v;
  typename AT::t_v_array d_omega;
  typename AT::t_f_array d_f;
  typename AT::t_f_array d_torque;
  typename AT::t_int_1d d_mask;
  typename AT::t_float_1d d_rmass;
  typename AT::t_float_1d d_radius;
  typename AT::tdual_float_2d k_history_one;
  typename AT::t_float_2d d_history_one;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix wall/gran requires atom style sphere

Self-explanatory.

E: Invalid fix wall/gran interaction style

UNDOCUMENTED

E: Cannot use wall in periodic dimension

Self-explanatory.

E: Cannot wiggle and shear fix wall/gran

Cannot specify both options at the same time.

E: Invalid wiggle direction for fix wall/gran

Self-explanatory.

E: Invalid shear direction for fix wall/gran

Self-explanatory.

E: Cannot wiggle or shear with fix wall/gran/region

UNDOCUMENTED

U: Fix wall/gran is incompatible with Pair style

Must use a granular pair style to define the parameters needed for
this fix.

*/
