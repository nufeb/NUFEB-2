#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build and install vtk-8.0.0 package.
#-------------------------------------------------------------------------------------------

cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD

vtk_dir="VTK-8.0.0"
vtk_file="VTK-8.0.0.tar.gz"
vtk_url="https://www.vtk.org/files/release/8.0/VTK-8.0.0.tar.gz"

download_source() {
    
  [ -d "$vtk_dir" ] &&  echo "$vtk_dir exists.  Not downloading..." && return 0

  echo "Downloading and extracting $vtk_dir"
  wget -O - $vtk_url | tar xvz
}
  	
download_source

cd $vtk_dir || exit 1
mkdir vtk-build
cd vtk-build || exit 1
mkdir vtk-8.0
intallpath=$PWD/vtk-8.0
echo $intallpath

cmake \
 -DBUILD_TESTING=OFF \
 -DCMAKE_INSTALL_PREFIX:PATH=$intallpath ../

make -j4
make install

version=`uname`
# set LD path according to different versions

if grep -q $intallpath ~/.bashrc; then
  echo -n
else
  echo "Writing path to .bashrc"
  if [ $version == "Linux" ]; then
    echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$intallpath/lib" >> ~/.bashrc
  elif [ $version == "Darwin" ]; then
    echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$intallpath/lib" >> ~/.bashrc
  fi
fi

exit 1
