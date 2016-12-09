/*
 * test_fix_division.cpp
 *
 *  Created on: 3 Dec 2016
 *      Author: bowen
 */

#include "test_fix_kinetics_monod.h"

using namespace std;

using namespace LAMMPS_NS;

namespace {

// The fixture for testing class Foo.
  class FixKineticsMonodTest : public ::testing::Test {
  protected:
    FixKineticsMonodTest() {
      kinetics = NULL;
      avec = NULL;

      char *argv[] = {"fix_kinetics_monod", "-screen", "none", "-log", "none", NULL};
      int narg = 5;
      // Create a LAMMPS Object
      lmp = new LAMMPS(narg,argv,MPI_COMM_WORLD);
    }

    virtual ~FixKineticsMonodTest(){
      lmp->destroy();
    }

    virtual void SetUp() {
      // Set atom_style
      lmp->input->one("atom_style bio");
      // Modify atom_style parameters
      lmp->input->one("atom_modify map array sort 1000 5.0e-7");
      lmp->input->one("read_data_bio inputs/monod.in");
      // Create group
      lmp->input->one("group HET type 1");
      // Define variables
      lmp->input->one("variable temp equal 298.15");
      lmp->input->one("variable gasTran equal 0.08205746");
      lmp->input->one("variable gvol equal 2e-11");
      lmp->input->one("variable EPSdens equal 30");

      avec = (AtomVecBio *) lmp->atom->style_match("bio");
    }

    virtual void TearDown() {
      lmp->input->one("clear");
    }

    // Objects declared here can be used by all tests in the test case for Foo.
    LAMMPS* lmp;
    AtomVecBio *avec;
    FixKinetics *kinetics;
  };

  /*
   *  Check for new mass/radius after non-diffusion growth
   */
 TEST_F(FixKineticsMonodTest, No_Diffusion_1step) {
   lmp->input->one("timestep 36000");
   lmp->input->one("fix k1 all kinetics 1 1 1 v_temp");
   lmp->input->one("fix ki1 all kinetics/monod 1 0 v_gasTran v_gvol v_EPSdens");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix >= 0) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   // run lammps
   lmp->input->one("run 1");

   EXPECT_NEAR(3.02e-17, lmp->atom->rmass[0], 1e-19);
   EXPECT_NEAR(3.63e-07, lmp->atom->radius[0], 1e-9);
   EXPECT_NEAR(7.68e-17, avec->outerMass[0], 1e-19);
   EXPECT_NEAR(8.70e-07, avec->outerRadius[0], 1e-9);
 }

 /*
  *  Check for new mass/radius after diffusion growth
  */
 TEST_F(FixKineticsMonodTest, Diffusion_1step_Consumption) {
   TearDown();
   SetUp();
   lmp->input->one("variable diffT equal 1e-3");
   lmp->input->one("variable tol equal 1e-6");

   lmp->input->one("timestep 36000");
   lmp->input->one("fix k1 all kinetics 1 1 1 v_temp");
   lmp->input->one("fix g1 all diffusion exp v_diffT v_tol pp pp nd");
   lmp->input->one("fix km1 all kinetics/monod 1 1 v_gasTran v_gvol v_EPSdens");

   int ifix = lmp->modify->find_fix("k1");
   if (ifix > -1) {
     kinetics = static_cast<FixKinetics *>(lmp->modify->fix[ifix]);
   }
   // run lammps
   lmp->input->one("run 1");

   EXPECT_NEAR(-9.01e-10, kinetics->nuR[1][0], 1e-12);
 }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
