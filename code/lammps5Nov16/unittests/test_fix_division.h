#ifndef FIX_DIVISION_TEST_H
#define FIX_DIVISION_TEST_H

#include "gtest/gtest.h"
#include <iostream>
#include <math.h>
#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "lammps.h"
#include "atom_vec.h"
#include "comm.h"
#include "input.h"
#include "integrate.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "mpi.h"

class FixDivisionTest{
protected:
  FixDivisionTest();
  virtual ~FixDivisionTest();
  virtual void SetUp();
  virtual void TearDown();
};

#endif
