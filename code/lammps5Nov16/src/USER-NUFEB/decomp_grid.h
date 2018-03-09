#ifndef LMP_DECOMP_GRID_H
#define LMP_DECOMP_GRID_H

#include "subgrid.h"

#include <vector>

namespace LAMMPS_NS {
template <class Derived>
class DecompGrid {
 public:
  void setup() {
    Derived *derived = static_cast<Derived *>(this);
    // communicate grid extent
    Grid<double, 3> grid = derived->get_grid();
    Subgrid<double, 3> subgrid = derived->get_subgrid();
    Box<int, 3> box = subgrid.get_box();
    std::vector<int> boxlo(3 * derived->comm->nprocs); 
    MPI_Allgather(&box.lower[0], 3, MPI_INT, boxlo.data(), 3, MPI_INT, derived->world);
    std::vector<int> boxhi(3 * derived->comm->nprocs);
    MPI_Allgather(&box.upper[0], 3, MPI_INT, boxhi.data(), 3, MPI_INT, derived->world);

    clear();
    recv_begin.resize(derived->comm->nprocs);
    send_begin.resize(derived->comm->nprocs);
    recv_end.resize(derived->comm->nprocs);
    send_end.resize(derived->comm->nprocs);
    requests.resize(derived->comm->nprocs);
    std::array<bool, 3> periodic = derived->get_periodic_boundary();
    // look for intersections
    int nrecv = 0;
    int nsend = 0;
    Box<int, 3> extbox = extend(box);
    Subgrid<double, 3> extgrid(grid, extbox);
    for (int p = 0; p < derived->comm->nprocs; p++) {
      recv_begin[p] = nrecv;
      send_begin[p] = nsend;
      Box<int, 3> other(&boxlo[3 * p], &boxhi[3 * p]);
      if (p != derived->comm->me) {
	// identify which cell we are going to send and receive
	setup_comm_cells(extgrid, box, other, nsend, nrecv);
	// check for periodic boundary conditions
	if (extbox.lower[0] < 0 && periodic[0]) {
	  setup_comm_cells(extgrid, box, translate(other, {-grid.get_dimensions()[0], 0, 0}), nsend, nrecv);
	}
	if (extbox.upper[0] > grid.get_dimensions()[0] && periodic[0]) {      
	  setup_comm_cells(extgrid, box, translate(other, {grid.get_dimensions()[0], 0, 0}), nsend, nrecv);
	}
	if (extbox.lower[1] < 0 && periodic[1]) {
	  setup_comm_cells(extgrid, box, translate(other, {0, -grid.get_dimensions()[1], 0}), nsend, nrecv);
	}
	if (extbox.upper[1] > grid.get_dimensions()[1] && periodic[1]) {      
	  setup_comm_cells(extgrid, box, translate(other, {0, grid.get_dimensions()[1], 0}), nsend, nrecv);
	}
	if (extbox.lower[2] < 0 && periodic[2]) {
	  setup_comm_cells(extgrid, box, translate(other, {0, 0, -grid.get_dimensions()[2]}), nsend, nrecv);
	}
	if (extbox.upper[2] > grid.get_dimensions()[2] && periodic[2]) {      
	  setup_comm_cells(extgrid, box, translate(other, {0, 0, grid.get_dimensions()[2]}), nsend, nrecv);
	}
      }
      recv_end[p] = nrecv;
      send_end[p] = nsend;
    }
    recv_buff.resize(derived->get_cell_data_size(nrecv));
    send_buff.resize(derived->get_cell_data_size(nsend));
  }

  void exchange() {
    Derived *derived = static_cast<Derived *>(this);
    // pack data to send buffer
    derived->pack_cells(send_cells.begin(), send_cells.end(), send_buff.begin());
    // send and recv grid data
    int nrequests = 0;
    int recv_offset = 0;
    int send_offset = 0;
    for (int p = 0; p < derived->comm->nprocs; p++) {
      if (p == derived->comm->me)
	continue;
      int nrecv = recv_end[p] - recv_begin[p];
      if (nrecv > 0) {
	int count = derived->get_cell_data_size(nrecv);
	MPI_Irecv(&recv_buff[recv_offset], count, MPI_DOUBLE, p, 0, derived->world, &requests[nrequests++]);
	recv_offset += count;
      }
      int nsend = send_end[p] - send_begin[p];
      if (nsend > 0) {
	int count = derived->get_cell_data_size(nsend);
	MPI_Isend(&send_buff[send_offset], count, MPI_DOUBLE, p, 0, derived->world, &requests[nrequests++]);
	send_offset += count;
      }
    }
    // wait for all MPI requests
    if (derived->comm->nprocs > 1) MPI_Waitall(nrequests, requests.data(), MPI_STATUS_IGNORE);
    // unpack data from recv buffer
    derived->unpack_cells(recv_cells.begin(), recv_cells.end(), recv_buff.begin());
  }

 private:
  void add_cells(const Subgrid<double, 3> &subgrid, const Box<int, 3> &box, std::vector<int> &cells)
  {
    for (int k = box.lower[2]; k < box.upper[2]; k++) {
      for (int j = box.lower[1]; j < box.upper[1]; j++) {
	for (int i = box.lower[0]; i < box.upper[0]; i++) {
	  cells.push_back(subgrid.get_linear_index({i, j, k}));
	}
      }
    }
  }

  bool check_intersection(const Box<int, 3> &g)
  {
    // check if the intersection is empty
    if (is_empty(g))
      return false;
    // check if the intersection is a corner
    int n[3];
    for (int i = 0; i < 3; i++)
      n[i] = g.upper[i] - g.lower[i];
    if ((n[0] > 1 && n[1] > 1) || (n[0] > 1 && n[2] > 1) || (n[1] > 1 && n[2] > 1))
      return true;
    return false;
  }

  int cell_count(const Box<int, 3> &box) {
    std::array<int, 3> s = size(box);
    return s[0] * s[1] * s[2];
  }

  void setup_comm_cells(const Subgrid<double, 3> &subgrid, const Box<int, 3> &box, const Box<int, 3> &other, int &nsend, int &nrecv) {
    // identify which cells we need to recv
    Box<int, 3> recvbox = intersect(extend(box), other);
    int n = cell_count(recvbox);
    if (check_intersection(recvbox)) {
      add_cells(subgrid, recvbox, recv_cells);
      nrecv += n;
    }
    // identify which cells we need to send
    Box<int, 3> sendbox = intersect(extend(other), box);
    n = cell_count(sendbox);
    if (check_intersection(sendbox)) {
      add_cells(subgrid, sendbox, send_cells);
      nsend += n;
    }
  }

  void clear() {
    recv_cells.clear();
    send_cells.clear();
    recv_begin.clear();
    recv_end.clear();
    send_begin.clear();
    send_end.clear();
    recv_buff.clear();
    send_buff.clear();
    requests.clear();
  }

  Grid<double, 3> grid;
  Grid<double, 3> subgrid;
  std::vector<int> recv_cells;
  std::vector<int> send_cells;
  std::vector<int> recv_begin;
  std::vector<int> recv_end;
  std::vector<int> send_begin;
  std::vector<int> send_end;
  std::vector<double> recv_buff;
  std::vector<double> send_buff;
  std::vector<MPI_Request> requests;
};
}

#endif // LMP_DECOMP_GRID_H
