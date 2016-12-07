/*
 * test_fix_division.cpp
 *
 *  Created on: 3 Dec 2016
 *      Author: bowen
 */

#include "test_fix_diffusion.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixDiffusionTest : public ::testing::Test {
  protected:
    FixDiffusionTest() {
      kinetics = NULL;
      diffusion = NULL;

      char *argv[] = {"fix_diffusion", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      // Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixDiffusionTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");
      lmp->input->one("read_data_bio inputs/kinetics.in");
      // Create group
      lmp->input->one("group HET type 1");
      // Define variables
      lmp->input->one("variable temp equal 298.15");
      lmp->input->one("variable diffT equal 1e-3");
      lmp->input->one("variable tol equal 1e-6");


    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    FixKinetics *kinetics;
    FixDiffusion *diffusion;
  };


 TEST_F(FixDiffusionTest, Explicit) {
   lmp->input->one("fix k1 all kinetics 5 5 5 v_temp");
   lmp->input->one("fix g1 all diffusion exp v_diffT v_tol pp pp pp");
   lmp->input->one("fix km1 all kinetics/monod 1 1 v_gasTran v_gvol v_EPSdens");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix >= 0) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   int dfix = lmp->modify->find_fix("g1");
   if (dfix >= 0) {
     diffusion = static_cast<FixDiffusion *>(lmp->modify->fix[dfix]);
   }
   // run lammps
   diffusion->diffusion(0);
 }


 TEST_F(FixDiffusionTest, Implicit) {
   TearDown();
   SetUp();
   lmp->input->one("fix k1 all kinetics 5 5 5 v_temp");
   lmp->input->one("fix g1 all diffusion imp v_diffT v_tol pp pp pp");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix >= 0) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   int dfix = lmp->modify->find_fix("g1");
   if (dfix >= 0) {
     diffusion = static_cast<FixDiffusion *>(lmp->modify->fix[dfix]);
   }
   // run lammps
   diffusion->diffusion(1);
 }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
