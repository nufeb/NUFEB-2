/*
 * test_fix_division.cpp
 *
 *  Created on: 3 Dec 2016
 *      Author: bowen
 */

#include "test_fix_kinetics_diffusion.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixDiffusionTest : public ::testing::Test {
  protected:
    FixDiffusionTest() {
      kinetics = NULL;

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
      // Set initial mass and growth rate to 0, no consumption
      lmp->input->one("read_data_bio inputs/diffusion.in");
      // Create group
      lmp->input->one("group HET type 1");
      // Define variables
      lmp->input->one("variable temp equal 298.15");
      lmp->input->one("variable diffT equal 1e-3");
      lmp->input->one("variable tol equal 1e-6");

      lmp->input->one("variable gvol equal 2e-11");
      lmp->input->one("variable gasTran equal 0.08205746");
      lmp->input->one("variable EPSdens equal 30");

      lmp->input->one("timestep 1");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    double getMinS(double **&nuS) {
      double result = 1;
      for (int i = 0; i < 125; i++) {
        if (nuS[1][i] < result)
          result = nuS[1][i];
      }
      return result;
    }

    double getMaxS(double **&nuS) {
      double result = 0;
      for (int i = 0; i < 125; i++) {
        if (nuS[1][i] > result)
          result = nuS[1][i];
      }
      return result;
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    FixKinetics *kinetics;
  };

  /*
   *  Check for nutrient concentration, use dirichlet BC for HighZ surface
   *  Top area has higher concentration than low area
   */
 TEST_F(FixDiffusionTest, Explicit_HighZ) {
   lmp->input->one("fix k1 all kinetics 5 5 5 v_temp");
   lmp->input->one("fix g1 all diffusion exp v_diffT v_tol pp pp nd");
   lmp->input->one("fix km1 all kinetics/monod 1 1 v_gasTran v_gvol v_EPSdens");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix >= 0) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   // run lammps
   lmp->input->one("run 1");
   double min = getMinS(kinetics->nuS);
   double max = getMaxS(kinetics->nuS);

   // Maximal value
   EXPECT_NEAR(0.029875, min, 1e-5);
   // Minimal value
   EXPECT_NEAR(0.029980, max, 1e-5);

   // Bottom area is supposed to be low S
   EXPECT_NEAR(0.029875, kinetics->nuS[1][0], 1e-5);
   EXPECT_NEAR(0.029875, kinetics->nuS[1][24], 1e-5);

   // Top area is supposed to be high S
   EXPECT_NEAR(0.029980, kinetics->nuS[1][100], 1e-5);
   EXPECT_NEAR(0.029980, kinetics->nuS[1][124], 1e-5);
 }


 /*
  *  Check for nutrient concentration, use dirichlet BC for LowZ surface
  *  Low area has higher concentration than high area
  */
 TEST_F(FixDiffusionTest, Explicit_LowZ) {
   TearDown();
   SetUp();
   lmp->input->one("fix k1 all kinetics 5 5 5 v_temp");
   lmp->input->one("fix g1 all diffusion exp v_diffT v_tol pp pp dn");
   lmp->input->one("fix km1 all kinetics/monod 1 1 v_gasTran v_gvol v_EPSdens");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix >= 0) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   // run lammps
   lmp->input->one("run 1");
   // Bottom surface is supposed to be high S
   EXPECT_NEAR(0.029980, kinetics->nuS[1][0], 1e-5);
   EXPECT_NEAR(0.029980, kinetics->nuS[1][24], 1e-5);

   // Top surface is supposed to be low S
   EXPECT_NEAR(0.029875, kinetics->nuS[1][100], 1e-5);
   EXPECT_NEAR(0.029875, kinetics->nuS[1][124], 1e-5);
 }


// TEST_F(FixDiffusionTest, Implicit) {
//   TearDown();
//   SetUp();
//   lmp->input->one("fix k1 all kinetics 5 5 5 v_temp");
//   lmp->input->one("fix g1 all diffusion imp v_diffT v_tol pp pp pp");
//   lmp->input->one("fix km1 all kinetics/monod 1 1 v_gasTran v_gvol v_EPSdens");
//
//   int ifix = lmp->modify->find_fix("k1");
//   if (ifix >= 0) {
//     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
//   }
//   // run lammps
//   lmp->input->one("run 1");
//   double min = getMinS(kinetics->nuS);
//   double max = getMaxS(kinetics->nuS);
//
//   // Maximal value
//   EXPECT_NEAR(0.029875, min, 1e-5);
//   // Minimal value
//   EXPECT_NEAR(0.029980, max, 1e-5);
//
//   // Bottom surface is supposed to be low S
//   EXPECT_NEAR(0.029875, kinetics->nuS[1][0], 1e-5);
//   EXPECT_NEAR(0.029875, kinetics->nuS[1][24], 1e-5);
//
//   // Top surface is supposed to be high S
//   EXPECT_NEAR(0.029980, kinetics->nuS[1][100], 1e-5);
//   EXPECT_NEAR(0.029980, kinetics->nuS[1][124], 1e-5);
// }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
