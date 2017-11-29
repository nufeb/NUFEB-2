/*
 * test_fix_division.cpp
 *
 *  Created on: 3 Dec 2016
 *      Author: bowen
 */

#include "test_fix_epsadh.h"
#include <string>
#include <sstream>

using namespace std;

using namespace LAMMPS_NS;

namespace {

  class FixEPSAdhTest : public ::testing::Test {
  protected:
    FixEPSAdhTest() {
      char *argv[] = {"fix_epsadh", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixEPSAdhTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");
      lmp->input->one("comm_modify  vel yes");
      // Read input data
      lmp->input->one("read_data_bio inputs/eps_adh.in");
      // Create groups
      lmp->input->one("group HET type 1");
      lmp->input->one("group EPS type 2");
      // Neighbor list
      lmp->input->one("neighbor 5.0e-7 bin");
      lmp->input->one("neigh_modify delay 0 one 5000");
      // Pair style
      lmp->input->one("pair_style  gran/hooke/history 0.000200 NULL 0.000050 NULL 0.0 1");
      lmp->input->one("pair_coeff  * *");

      lmp->input->one("timestep 1");
      // Define variables
      lmp->input->one("variable ke equal 5e+10");
      lmp->input->one("fix j1 all epsadh 1 v_ke 1");

      lmp->input->one("run 1");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS *lmp;
  };

  /*
   *  Check for adhesion force in z-axis
   */
 TEST_F(FixEPSAdhTest, EPS_AxisZ) {

   EXPECT_NEAR(0, lmp->atom->f[0][0], 1e-11);
   EXPECT_NEAR(0, lmp->atom->f[0][1], 1e-11);
   EXPECT_NEAR(-1.525323e-09, lmp->atom->f[0][2], 1e-11);

  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
