#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build and install hdf5-1.10.5 package.
#-------------------------------------------------------------------------------------------


set -euo pipefail

cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD

dir_name="hdf5-1.10.5"
hdf5_file="hdf5-1.10.5.tar.gz"
hdf_url="https://www.hdfgroup.org/package/hdf5-1-10-5-tar-gz/?wpdmdl=13571&refresh=63eb7b4c5b2e71676376908"

download_source() {
    
  [ -d "$dir_name" ] &&  echo "$dir_name exists.  Not downloading..." && return 0

  echo "Downloading $dir_name"
  wget -O $hdf5_file $hdf_url 
  
  echo "Extracting $hdf5_file..."
  tar xzf $hdf5_file
  
  echo "Deleting $hdf5_file..."
  rm $hdf5_file
}
  	
download_source

cd hdf5-1.10.5 || exit 1

./configure --enable-parallel --enable-shared 

echo "Building $dir_name..."
make
make install

version=`uname`

# set LD path according to different versions
intallpath=$currentDir/$dir_name/hdf5

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

echo 'Done'

exit 1
