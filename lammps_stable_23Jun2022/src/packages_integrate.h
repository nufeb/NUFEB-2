#undef PACKAGE
#define PACKAGE "INTEL"
#include "INTEL/verlet_lrt_intel.h"
#undef PACKAGE
#define PACKAGE "KOKKOS"
#include "KOKKOS/nufeb_run_kokkos.h"
#undef PACKAGE
#define PACKAGE "KOKKOS"
#include "KOKKOS/verlet_kokkos.h"
#undef PACKAGE
#define PACKAGE "OPENMP"
#include "OPENMP/respa_omp.h"
#undef PACKAGE
#define PACKAGE "REPLICA"
#include "REPLICA/verlet_split.h"
#undef PACKAGE
#define PACKAGE "USER-NUFEB"
#include "USER-NUFEB/nufeb_run.h"
