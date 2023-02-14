#!/bin/bash

set -euo pipefail

cd ${0%/*} || exit 1 # Run from this directory

echo "Uninstalling NUFEB.."
currentDir=$PWD

cd $currentDir/lammps_stable_29Oct2020/src || exit 1
make no-all
make no-kokkos
make clean-all
