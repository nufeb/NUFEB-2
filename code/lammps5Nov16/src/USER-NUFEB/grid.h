#ifndef SRC_GRID_H
#define SRC_GRID_H

#include <array>

namespace LAMMPS_NS {

template <typename T, std::size_t N>
struct Grid {
  Grid() {
    for (int i = 0; i < N; i++) {
      lower[i] = 0;
      upper[i] = 0; 
    }    
  }
  Grid(int *l, int *u) {
    for (int i = 0; i < N; i++) {
      lower[i] = l[i];
      upper[i] = u[i]; 
    }
  }

  T lower[N];
  T upper[N];	
};

template <typename T, std::size_t N>
inline Grid<T, N> intersect(const Grid<T, N> &g0, const Grid<T, N> &g1) {
  Grid<T, N> ret;
  for (int i = 0; i < N; i++) {
    ret.lower[i] = g0.lower[i] > g1.lower[i] ? g0.lower[i] : g1.lower[i];
    ret.upper[i] = g0.upper[i] < g1.upper[i] ? g0.upper[i] : g1.upper[i];
  }
  return ret;
}

template <typename T, std::size_t N>
inline Grid<T, N> extend(const Grid<T, N> &g, T n = T(1))
{
  Grid<T, N> ret;
  for (int i = 0; i < N; i++) {
    ret.lower[i] = g.lower[i] - n;
    ret.upper[i] = g.upper[i] + n;
  }
  return ret;
}

template <typename T, std::size_t N>
inline bool is_empty(const Grid<T, N> &g)
{
  for (int i = 0; i < 3; i++) {
    if (g.lower[i] >= g.upper[i])
      return true;
  }
  return false;
}

template <typename T, std::size_t N>
inline T get_linear_index(const Grid<T, N> &g, std::array<T, N> cell)
{
  T n = T(1);
  T result = T(0);
  for (int i = 0; i < N; i++) {
    result += (cell[i] - g.lower[i]) * n;
    n *= g.upper[i] - g.lower[i];
  }
  return result;
}

template <typename T, std::size_t N>
inline int cell_count(const Grid<T, N> &g)
{
  T result = T(1);
  for(int i = 0; i < N; i++) {
    result *= g.upper[i] - g.lower[i];
  }
  return result;
}

template <typename T, std::size_t N>
inline Grid<T, N> translate(const Grid<T, N> &g, std::array<T, N> d) {
  Grid<T, N> ret(g);
  for (int i = 0; i < N; i++) {
    ret.lower[i] += d[i];
    ret.upper[i] += d[i];
  }
  return ret;
}

} // namespace LAMMPS_NS

#endif // SRC_GRID_H

