#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build and install openfoam-4.0
# The installation directory will be $HOME/OpenFOAM
#-------------------------------------------------------------------------------------------

cd ${0%/*} || exit 1 # Run from this directory

current_dir=$PWD
of4_dir="OpenFOAM"
inst_dir="$of4_dir/OpenFOAM-4.0"
third_dir="ThirdParty-4.0"

download_source() {
  [ -d "$of4_dir" ] &&  echo "$of4_dir exists.  Not downloading..." && return 0

  mkdir $of4_dir

  echo "Downloading and extracting $of4_dir"
  wget -O - http://dl.openfoam.org/source/4-0 | tar xvz
  wget -O - http://dl.openfoam.org/third-party/4-0 | tar xvz

  mv OpenFOAM-4.x-version-4.0 $inst_dir
  mv ThirdParty-4.x-version-4.0 $of4_dir/ThirdParty-4.0
}

download_source

cd $inst_dir || exit 1
echo "Building $of4_dir..."

arch=`uname -m`
echo $arch
if [ "$arch" = "i686" ]; then
  . $current_dir/$inst_dir/etc/bashrc WM_ARCH_OPTION=32 FOAMY_HEX_MESH=yes
elif [ "$arch" = "x86_64" ]; then
  . $current_dir/$inst_dir/etc/bashrc WM_LABEL_SIZE=64 FOAMY_HEX_MESH=yes
fi

./Allwmake -j 4
