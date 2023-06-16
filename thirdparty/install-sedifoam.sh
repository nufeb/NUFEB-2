#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build and install sediFoam
#-------------------------------------------------------------------------------------------

cd ${0%/*} || exit 1 # Run from this directory

sedi_dir="sediFoam"
sedi_url="https://github.com/xiaoh/sediFoam/tarball/master"

download_source() {

  [ -d "$sedi_dir" ] &&  echo "$sedi_dir exists.  Not downloading..." && return 0

  echo "Downloading and extracting $sedi_dir"
  wget -O - $sedi_url | tar xvz

  mv xiaoh-sediFoam* sediFoam
}

download_source