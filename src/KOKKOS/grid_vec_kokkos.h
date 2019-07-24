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

#ifndef LMP_GRID_VEC_KOKKOS_H
#define LMP_GRID_VEC_KOKKOS_H

#include "grid_vec.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {
class GridVecKokkos : public GridVec {
 public:
  GridVecKokkos(class LAMMPS *);
  virtual ~GridVecKokkos() {}
  virtual void setup();
  virtual void sync(ExecutionSpace, unsigned int) = 0;
  virtual void modified(ExecutionSpace, unsigned int) = 0;
  virtual void sync_overlapping_device(ExecutionSpace, unsigned int) = 0;
  virtual int pack_comm_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &) = 0;
  virtual void unpack_comm_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &) = 0;
  virtual int pack_exchange_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &) = 0;
  virtual void unpack_exchange_kokkos(int, int, const DAT::tdual_int_1d &, const DAT::tdual_xfloat_1d &) = 0;

 protected:
  size_t buffer_size;
  void* buffer;

  #ifdef KOKKOS_ENABLE_CUDA
  template<class ViewType>
  Kokkos::View<typename ViewType::data_type,
               typename ViewType::array_layout,
               Kokkos::CudaHostPinnedSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged> >
  create_async_copy(const ViewType& src) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same<typename ViewType::execution_space,LMPDeviceType>::value,
                   Kokkos::CudaHostPinnedSpace,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer = Kokkos::kokkos_malloc<Kokkos::CudaHostPinnedSpace>(src.span());
       buffer_size = src.span();
    } else if (buffer_size < src.span()) {
       buffer = Kokkos::kokkos_realloc<Kokkos::CudaHostPinnedSpace>(buffer,src.span());
       buffer_size = src.span();
    }
    return mirror_type( buffer ,
                             src.extent(0) ,
                             src.extent(1) ,
                             src.extent(2) ,
                             src.extent(3) ,
                             src.extent(4) ,
                             src.extent(5) ,
                             src.extent(6) ,
                             src.extent(7) );
  }

  template<class ViewType>
  void perform_async_copy(ViewType& src, unsigned int space) {
    typedef Kokkos::View<typename ViewType::data_type,
                 typename ViewType::array_layout,
                 typename std::conditional<
                   std::is_same<typename ViewType::execution_space,LMPDeviceType>::value,
                   Kokkos::CudaHostPinnedSpace,typename ViewType::memory_space>::type,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged> > mirror_type;
    if (buffer_size == 0) {
       buffer = Kokkos::kokkos_malloc<Kokkos::CudaHostPinnedSpace>(src.span()*sizeof(typename ViewType::value_type));
       buffer_size = src.span();
    } else if (buffer_size < src.span()) {
       buffer = Kokkos::kokkos_realloc<Kokkos::CudaHostPinnedSpace>(buffer,src.span()*sizeof(typename ViewType::value_type));
       buffer_size = src.span();
    }
    mirror_type tmp_view( (typename ViewType::value_type*)buffer ,
                             src.extent(0) ,
                             src.extent(1) ,
                             src.extent(2) ,
                             src.extent(3) ,
                             src.extent(4) ,
                             src.extent(5) ,
                             src.extent(6) ,
                             src.extent(7) );
    if(space == Device) {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.h_view),
      Kokkos::deep_copy(LMPHostType(),src.d_view,tmp_view);
      src.clear_sync_state();
    } else {
      Kokkos::deep_copy(LMPHostType(),tmp_view,src.d_view),
      Kokkos::deep_copy(LMPHostType(),src.h_view,tmp_view);
      src.clear_sync_state();
    }
  }
  #else
  template<class ViewType>
  void perform_async_copy(ViewType& src, unsigned int space) {
    if(space == Device)
      src.template sync<LMPDeviceType>();
    else
      src.template sync<LMPHostType>();
  }
  #endif
};
 
}

#endif
