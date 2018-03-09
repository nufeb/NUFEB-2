#ifndef LMP_GRID_H
#define LMP_GRID_H

#include "box.h"

#include <type_traits>

namespace LAMMPS_NS {
template <typename T, std::size_t N, typename IndexType = int>
class Grid {
 public:
  Grid() : origin({0}), dimensions({0}), cell_size({1}) {}
  Grid(const Grid &other) = default;
  Grid(const std::array<T, N> &origin,
    const std::array<IndexType, N> &dimensions,
    const std::array<T, N> &cell_size) :
    origin(origin), dimensions(dimensions), cell_size(cell_size) {}
  Grid(const Box<T, N> &b, const std::array<IndexType, N> &dimensions) :
    origin(b.lower), dimensions(dimensions) {
    for (int i = 0; i < N; i++) {
      cell_size[i] = (b.upper[i] - b.lower[i]) / dimensions[i];
    }
  }
  template<typename U>
 Grid(const Box<U, N> &b, typename std::enable_if<std::is_integral<U>::value && std::is_same<T, U>::value>::type* = 0) : Grid(b.lower, size(b), {1}) {}

  const std::array<T, N> &get_origin() const {
    return origin;
  }

  const std::array<IndexType, N> &get_dimensions() const {
    return dimensions;
  }

  const std::array<T, N> &get_cell_size() const {
    return cell_size;
  }

  IndexType get_linear_index(const std::array<IndexType, N> &cell) const {
    IndexType n = T(1);
    IndexType result = T(0);
    for (int i = 0; i < N; i++) {
      result += cell[i] * n;
      n *= dimensions[i];
    }
    return result;
  }

  IndexType cell_count() const {
    IndexType result = T(1);
    for(int i = 0; i < N; i++) {
      result *= dimensions[i];
    }
    return result;
  }

  std::array<IndexType, N> get_index(const std::array<T, N> &x) const {
    std::array<IndexType, N> result;
    for (int i = 0; i < N; i++) {
      result[i] = (x[i] - origin[i]) / cell_size[i];
    }
    return result;
  }

 private:
  std::array<T, N> origin;
  std::array<IndexType, N> dimensions;
  std::array<T, N> cell_size;
  Box<T, N> box;
};
} // namespace LAMMPS_NS
#endif // LMP_GRID_H
