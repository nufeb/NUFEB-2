/*
 * test_fix_division.cpp
 *
 *  Created on: 3 Dec 2016
 *      Author: bowen
 */

#include "test_fix_division.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixDivisionTest : public ::testing::Test {
  protected:
    FixDivisionTest() {
      char *argv[] = {"fix_division", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      //Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");

      lmp->input->one("comm_modify vel yes");
    }

    virtual ~FixDivisionTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      //Initialize simulation domain, type, atom and nutrients
      // Read input data
      lmp->input->one("read_data_bio one_atom.in");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
  };

  // Check for case with to high a value of bcoeff_narg
 TEST_F(FixDivisionTest, mytest) {
    ASSERT_TRUE(lmp->atom->ntypes == 1);
 }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
