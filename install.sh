#!/bin/bash

set -euo pipefail

cd ${0%/*} || exit 1 # Run from this directory

echo "Installing NUFEB.."
rootDir=$PWD

#### Copy package and lib files to LAMMPS directory #####
echo "Copying packages to LAMMPS.."
cp -rpf $rootDir/src/* $rootDir/lammps_stable_29Oct2020/src/
cp -rpf $rootDir/lib/* $rootDir/lammps_stable_29Oct2020/lib/

echo "Configuring Makefile.lammps.."

cd $rootDir/lammps_stable_29Oct2020/lib/nufeb || exit 1 
cp Makefile.lammps_core Makefile.lammps

declare -i vtk_hdf=0

for var in "$@"
do 
    if [ $var == "--enable-vtk" ] ; then
       cp Makefile.lammps_vtk8.0 Makefile.lammps
       cp ../vtk/Makefile.lammps_vtk8.0 ../vtk/Makefile.lammps
       vtk_hdf=$((vtk_hdf+1))
    elif [ $var == "--enable-hdf5" ]; then
       cp Makefile.lammps_hdf5 Makefile.lammps
       vtk_hdf=$((vtk_hdf+1))
    elif [ $var == "--enable-misc" ]; then continue
    elif [ $var == "--enable-plasmid" ]; then continue
    elif [ $var == "--gpu" ]; then continue
    elif [ $var == "--static" ]; then continue
    elif [ $var == "--shared" ]; then continue
    elif [ $var == "--serial" ]; then continue
    else
       echo "Unknown parameter: $var"
       exit 1
    fi
done

if [ $vtk_hdf = 2 ]; then
echo $PWD
    cp Makefile.lammps_hdf5_vtk8.0 Makefile.lammps
fi

#### Build LAMMPS with NUFEB and user defined packages#####
echo "Installing required packages.."

cd $rootDir/lammps_stable_29Oct2020/src || exit 1
make yes-user-nufeb
make yes-granular

for var in "$@"
do 
    if [ $var == "--enable-vtk" ]; then
	make yes-user-vtk
    elif [ $var == "--enable-hdf5" ]; then
	make yes-user-hdf5
    elif [ $var == "--enable-misc" ]; then
	make yes-misc
    elif [ $var == "--enable-plasmid" ]; then
	make yes-user-plasmid
    elif [ $var == "--gpu" ]; then
	make yes-kokkos
    fi
done

#### Write path to .bashrc#####
#if grep -q  "export PATH=\$PATH:$rootDir" ~/.bashrc; then
#   echo -n
#else
#   echo "Writing NUFEB root path to .bashrc"
#   echo "export PATH=\$PATH:$rootDir" >> ~/.bashrc
#fi


echo "Building NUFEB.."
for var in "$@"
do 
    if [ $var == "--serial" ]; then
	cd STUBS || exit 1
        make
        cd .. || exit 1;
        make -j4 serial
        mv lmp_serial $rootDir/nufeb_serial
        exit 1
    elif [ $var == "--static" ]; then   
        make -j4 mpi mode=lib
        exit 1
    elif [ $var == "--shared" ]; then   
        make -j4 mpi mode=shlib
        exit 1
    elif [ $var == "--gpu" ]; then   
        make -j4 kokkos_cuda_mpi
        mv lmp_kokkos_cuda_mpi $rootDir/nufeb_gpu
        exit 1 
   fi
done

make -j4 mpi
mv lmp_mpi $rootDir/nufeb_mpi
exit 1

