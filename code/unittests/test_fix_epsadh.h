#ifndef FIX_EPSADH_TEST_H
#define FIX_EPSADH_TEST_H

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
#include "mpi.h"
#include <string>
#include <sstream>

class FixEPSAdhTest{
protected:
  FixEPSAdhTest();
  virtual ~FixEPSAdhTest();
  virtual void SetUp();
  virtual void TearDown();
};

#endif
