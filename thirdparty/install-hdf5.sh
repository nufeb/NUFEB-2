#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build and install hdf5-1.10.5 package.
#-------------------------------------------------------------------------------------------

cd ${0%/*} || exit 1 # Run from this directory

currentDir=$PWD

hdf5_dir="hdf5-1.10.5"
hdf5_url="https://www.hdfgroup.org/package/hdf5-1-10-5-tar-gz/?wpdmdl=13571&refresh=63eb7b4c5b2e71676376908"

download_source() {
    
  [ -d "$hdf5_dir" ] &&  echo "$hdf5_dir exists.  Not downloading..." && return 0

  echo "Downloading and extracting $hdf5_dir"
  wget -O - $hdf5_url | tar xvz
}
  	
download_source

cd hdf5-1.10.5 || exit 1

./configure --enable-parallel --enable-shared 

echo "Building $hdf5_dir..."
make
make install

version=`uname`

# set LD path according to different versions
intallpath=$currentDir/$hdf5_dir/hdf5

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
