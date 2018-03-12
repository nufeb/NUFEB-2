#ifndef LMP_SUBGRID_H
#define LMP_SUBGRID_H

#include "grid.h"

namespace LAMMPS_NS {
template <typename T, std::size_t N, typename IndexType = int>
class Subgrid {
 public:
  Subgrid() = default;
  Subgrid(const Subgrid &other) = default;
  Subgrid(const Grid<T, N, IndexType> &grid, const std::array<IndexType, N> &origin, const std::array<IndexType, N> &dimensions) :
    grid(grid), origin(origin), dimensions(dimensions) {}
  Subgrid(const Grid<T, N, IndexType> &grid, const Box<IndexType, N> & box) :
    Subgrid(grid, box.lower, size(box)) {}
  Subgrid(const Grid<T, N, IndexType> &grid, const Box<T, N> & box) :
    grid(grid) {
    for (int i = 0; i < N; i++) {
      origin[i] = static_cast<IndexType>(box.lower[i] / grid.get_cell_size()[i]);
      dimensions[i] = static_cast<IndexType>(box.upper[i] / grid.get_cell_size()[i]) - origin[i];
    }
  }

  const Grid<T, N, IndexType> &get_grid() const {
    return grid;
  }

  const std::array<IndexType, N> &get_origin() const {
    return origin;
  }
 
  const std::array<IndexType, N> &get_dimensions() const {
    return dimensions;
  }

  Box<IndexType, N> get_box() const {
    Box<IndexType, N> result;
    result.lower = origin;
    for (int i = 0; i < N; i++) {
      result.upper[i] = origin[i] + dimensions[i];
    }
    return result;
  }

  IndexType get_linear_index(const std::array<IndexType, N> &cell) const {
    IndexType n = T(1);
    IndexType result = T(0);
    for (int i = 0; i < N; i++) {
      result += (cell[i] - origin[i]) * n;
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

  IndexType get_index(const std::array<T, N> &x) const {
    IndexType n = T(1);
    IndexType result = T(0);
    for (int i = 0; i < N; i++) {
      result += (static_cast<IndexType>(x[i] / grid.get_cell_size()[i]) - origin[i]) * n;
      n *= dimensions[i];
    }
    return result;
  }

  bool is_inside(const std::array<T, N> &x) const {
    bool result = true;
    for (int i = 0; i < N; i++) {
      result &= x[i] >= origin[i] * grid.get_cell_size()[i] && x[i] < (origin[i] + dimensions[i]) * grid.get_cell_size()[i];  
    }
    return result;
  }

 private:
  Grid<T, N, IndexType> grid;
  std::array<IndexType, N> origin;
  std::array<IndexType, N> dimensions;
};
} // namespace LAMMPS_NS
#endif // LMP_SUBGRID_H
