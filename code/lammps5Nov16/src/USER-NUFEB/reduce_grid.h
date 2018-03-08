#ifndef LMP_REDUCE_GRID_H
#define LMP_REDUCE_GRID_H

#include "grid.h"

#include <vector>

namespace LAMMPS_NS {
template <class Derived>
class ReduceGrid {
 public:
  void setup() {
    Derived *derived = static_cast<Derived *>(this);
    // communicate grid extent
    Box<int, 2> box = derived->get_subgrid();
    Grid<int, 2> basegrid(box);
    std::vector<int> boxlo(2 * derived->comm->nprocs); 
    MPI_Allgather(&box.lower[0], 2, MPI_INT, boxlo.data(), 2, MPI_INT, derived->world);
    std::vector<int> boxhi(2 * derived->comm->nprocs);
    MPI_Allgather(&box.upper[0], 2, MPI_INT, boxhi.data(), 2, MPI_INT, derived->world);

    clear();
    recv_begin.resize(derived->comm->nprocs);
    send_begin.resize(derived->comm->nprocs);
    recv_end.resize(derived->comm->nprocs);
    send_end.resize(derived->comm->nprocs);
    requests.resize(derived->comm->nprocs);
    // look for intersections
    int nrecv = 0;
    int nsend = 0;
    for (int p = 0; p < derived->comm->nprocs; p++) {
      recv_begin[p] = nrecv;
      send_begin[p] = nsend;
      Box<int, 2> other(&boxlo[2 * p], &boxhi[2 * p]);
      if (p != derived->comm->me) {
	// identify which cells we need to recv if we are the bottom most proc
	Box<int, 2> intersection = intersect(box, other);
	int n = cell_count(intersection);
	if (n > 0 && derived->is_bottom_most()) {
	  add_cells(basegrid, intersection, recv_cells);
	  nrecv += n;
	}
	// identify which cells we need to send if we are not the bottom most proc
	if (n > 0 && !derived->is_bottom_most()) {
	  add_cells(basegrid, intersection, send_cells);
	  nsend += n;
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
    }
    for (int i = 0; i < nrequests; i++) {
      MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
      derived->unpack_cells_reduce(recv_cells.begin(), recv_cells.end(), recv_buff.begin(), [](auto v0, auto v1) { return MAX(v0, v1); });
    }
    // pack data to send buffer
    derived->pack_cells(send_cells.begin(), send_cells.end(), send_buff.begin());
    nrequests = 0;
    for (int p = 0; p < derived->comm->nprocs; p++) {
      if (p == derived->comm->me)
	continue;
      int nsend = send_end[p] - send_begin[p];
      if (nsend > 0) {
	int count = derived->get_cell_data_size(nsend);
	MPI_Isend(&send_buff[send_offset], count, MPI_DOUBLE, p, 0, derived->world, &requests[nrequests++]);
	send_offset += count;
      }
    }
    // wait for all MPI requests
    if (derived->comm->nprocs > 1) MPI_Waitall(nrequests, requests.data(), MPI_STATUS_IGNORE);
  }

 private:
  void add_cells(const Grid<int, 2> &basegrid, const Box<int, 2> &box, std::vector<int> &cells)
  {
    for (int j = box.lower[1]; j < box.upper[1]; j++) {
      for (int i = box.lower[0]; i < box.upper[0]; i++) {
	cells.push_back(basegrid.get_linear_index({i - basegrid.get_origin()[0], j - basegrid.get_origin()[1]}));
      }
    }
  }
  
  bool is_bottom_most() {
    Derived *derived = static_cast<Derived *>(this);
    return derived->domain->sublo[2] == derived->domain->boxlo[2];
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

#endif // LMP_REDUCE_GRID_H
