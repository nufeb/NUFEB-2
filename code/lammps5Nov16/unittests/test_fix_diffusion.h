#ifndef FIX_DIFFUSION_TEST_H
#define FIX_DIFFUSION_TEST_H

#include "gtest/gtest.h"
#include <iostream>
#include <math.h>
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "lammps.h"
#include "atom_vec.h"
#include "fix_kinetics.h"
#include "bio.h"
#include "comm.h"
#include "input.h"
#include "mpi.h"
#include "modify.h"
#include "update.h"
#include "library.h"

class FixDiffusionTest{
protected:
  FixDiffusionTest();
  virtual ~FixDiffusionTest();
  virtual void SetUp();
  virtual void TearDown();
};

#endif
