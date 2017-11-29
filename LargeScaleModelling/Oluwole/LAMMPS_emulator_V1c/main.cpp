// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// Creating a lattice plot from RInside
// cf http://stackoverflow.com/questions/24378223/saving-lattice-plots-with-rinside-and-rcpp/
//
// Copyright (C) 2014  Dirk Eddelbuettel and GPL'ed 

#include <RInside.h>               // for the embedded R via RInside
#include <unistd.h>

int main(int argc, char *argv[]) {

    // create an embedded R instance
    RInside R(argc, argv);               

    // evaluate an R expression using source command
    std::string cmd = "source('emu.R',echo=FALSE)";
  // Rcpp::NumericMatrix tmpfile = R.parseEval(cmd);
 //   std::string tmpfile = 
    R.parseEval(cmd);
 //   std::cout << tmpfile << std::endl;
     exit(0);
}
