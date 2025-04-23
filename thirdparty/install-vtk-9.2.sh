#!/bin/bash

#-------------------------------------------------------------------------------------------
# This script will download, configure, build, and install VTK 9.2.6 package.
#-------------------------------------------------------------------------------------------

cd "${0%/*}" || exit 1 # Run from this script's directory

currentDir=$PWD

vtk_dir="VTK-9.2.6"
vtk_repo="https://gitlab.kitware.com/vtk/vtk.git"
vtk_branch="v9.2.6"

# Download source if not already present
download_source() {
  if [ -d "$vtk_dir" ]; then
    echo "$vtk_dir already exists. Skipping download..."
    return 0
  fi

  echo "Cloning VTK from GitLab..."
  git clone --branch $vtk_branch --depth 1 $vtk_repo $vtk_dir
}

download_source

cd "$vtk_dir" || exit 1
mkdir -p build install
cd build || exit 1

install_path="$PWD/../install"
echo "Installing to $install_path"

cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$install_path" \
  -DBUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF \
  -DVTK_GROUP_ENABLE_Qt=NO \
  -DBUILD_SHARED_LIBS=ON

make -j$(nproc)
make install

# Update environment variables
if grep -q "$install_path" ~/.bashrc; then
  echo "Install path already in .bashrc"
else
  echo "Writing VTK install path to .bashrc"
  echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$install_path/lib" >> ~/.bashrc
  echo "export PATH=\$PATH:$install_path/bin" >> ~/.bashrc
fi

echo "VTK 9.2.6 installation complete!"
exit 0

