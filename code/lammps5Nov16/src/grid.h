#ifndef SRC_GRID_H
#define SRC_GRID_H

namespace LAMMPS_NS {

struct Grid
{
  Grid() {
    for (int i = 0; i < 3; i++) {
      lower[i] = 0;
      upper[i] = 0; 
    }    
  }
  Grid(int *l, int *u) {
  	for (int i = 0; i < 3; i++) {
      lower[i] = l[i];
      upper[i] = u[i]; 
  	}
  }

  int lower[3];
  int upper[3];	
};

inline Grid intersect(const Grid &g0, const Grid &g1)
{
  Grid ret;
  for (int i = 0; i < 3; i++) {
    ret.lower[i] = g0.lower[i] > g1.lower[i] ? g0.lower[i] : g1.lower[i];
    ret.upper[i] = g0.upper[i] < g1.upper[i] ? g0.upper[i] : g1.upper[i];
  }
  return ret;
}

inline Grid extend(const Grid &g)
{
	Grid ret;
 	for (int i = 0; i < 3; i++) {
 	  ret.lower[i] = g.lower[i] - 1;
    ret.upper[i] = g.upper[i] + 1;
  }
  return ret;
}

inline bool is_empty(const Grid &g)
{
  for (int i = 0; i < 3; i++) {
    if (g.lower[i] >= g.upper[i])
      return true;
  }
  return false;
}

inline int get_linear_index(const Grid &g, int i, int j, int k)
{
  int nx = g.upper[0] - g.lower[0];
  int ny = g.upper[1] - g.lower[1];
  return i - g.lower[0] + (j - g.lower[1]) * nx + (k - g.lower[2]) * nx * ny;
}

inline int cell_count(const Grid &g)
{
  int nx = g.upper[0] - g.lower[0];
  int ny = g.upper[1] - g.lower[1];
  int nz = g.upper[2] - g.lower[2];
  return nx * ny * nz; 
}

inline Grid translate(const Grid &g, int x, int y, int z) {
  Grid ret(g);
  ret.lower[0] += x;
  ret.lower[1] += y;
  ret.lower[2] += z;
  ret.upper[0] += x;
  ret.upper[1] += y;
  ret.upper[2] += z;
  return ret;
}

} // namespace LAMMPS_NS

#endif // SRC_GRID_H

