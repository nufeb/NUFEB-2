#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build and install sediFoam
#-------------------------------------------------------------------------------------------

cd ${0%/*} || exit 1 # Run from this directory

# Read the information of current directory.
# And collect information of the installation of LAMMPS from user.
echo "Installing sediFoam (for mac/linux).."

echo "*******************************************"
echo "select the system you are running, then press enter"
echo "  1) Ubuntu14.x - Ubuntu16.x"
echo "  2) Ubuntu17.x - Ubuntu18.x" 
echo "  3) Centos"
echo "  4) Mac" 
echo "*******************************************"
read n

case $n in
  1) echo "You chose 1) Ubuntu14.x - Ubuntu16.x";;
  2) echo "You chose 2) Ubuntu17.x - Ubuntu18.x";;
  3) echo "You chose 3) Centos";;
  4) echo "You chose 4) Mac";;
  *) echo "Unknown option"; exit;;
esac

sedi_dir="sediFoam"
sedi_url="https://github.com/nufeb/sediFoam/tarball/master"

download_source() {

  [ -d "$sedi_dir" ] &&  echo "$sedi_dir exists.  Not downloading..." && return 0

  echo "Downloading and extracting $sedi_dir"
  wget -O - $sedi_url | tar xvz

  mv nufeb-sediFoam* sediFoam
}

download_source

sedifoam_dir=$(readlink -f sediFoam)
nufeb_dir=$(readlink -f ..)
lammps_dir=$nufeb_dir/lammps_stable_23Jun2022
lammps_src=$lammps_dir/src

echo "Directory of LAMMPS is: " $lammps_src

echo "Copying packages to LAMMPS.."
cp -rf $sedifoam_dir/interfaceToLammps/* $lammps_src/

# Make STUBS 
cd $lammps_src/STUBS || exit 1
make
cd $lammps_src || exit 1

# Make packages
make yes-GRANULAR
make yes-NUFEB
make yes-COLLOID
make yes-SEDIFOAM

# Build NUFEB dynamic library
make -j4 shanghaimac mode=shlib
cd $FOAM_LIBBIN || exit 1
    
# Use different options according to different versions
if [ $n == 4 ];then 
    ln -sf $lammps_dir/src/liblammps_shanghaimac.so .
    cd $sedifoam_dir/lammpsFoam || exit 1
    touch Make/options
    echo "LAMMPS_DIR ="$lammps_src > Make/options
    cat Make/options-mac-openmpi >> Make/options
else
    ln -sf $lammps_dir/src/liblammps_shanghailinux.so .
    cd $sedifoam_dir/lammpsFoam || exit 1
    touch Make/options
    echo "LAMMPS_DIR ="$lammps_src > Make/options

    case $n in
	1) cat Make/options-ubuntu16-openmpi >> Make/options ;;
	2) cat Make/options-ubuntu18-openmpi >> Make/options ;;
	3) cat Make/options-linux-openmpi >> Make/options ;;
    esac

fi 

wmake libso dragModels 
wmake libso chPressureGrad 
wmake libso lammpsFoamTurbulenceModels
wmake 
