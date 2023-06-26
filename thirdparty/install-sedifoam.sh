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

sedifoam_dir=$PWD/sediFoam
cd .. || exit 1
nufeb_dir=$PWD
lammps_dir=$PWD/lammps_stable_23Jun2022
lammps_src=$lammps_dir/src

echo "Directory of LAMMPS is: " $lammpsDir

echo "Copying packages to LAMMPS.."
echo $lammps_src/
cp -rf $sedifoam_dir/interfaceToLammps/* $lammps_src/
echo $sedifoam_dir/interfaceToLammps/
cp -rf $nufeb_dir/src/* $lammps_src/
cp -rf $nufeb_dir/lib/* $lammps_src/lib/

# Make STUBS
cd $lammps_src/STUBS || exit 1
make
cd $lammps_src || exit 1

# Make packages
make yes-GRANULAR
make yes-NUFEB
make yes-COLLOID
make yes-SEDIFOAM

#make -j4 shanghailinux mode=shared
cd $FOAM_USER_LIBBIN || exit 1
ln -sf $lammps_dir/src/liblammps_shanghailinux.so .
cd $sedifoam_dir/lammpsFoam || exit 1
touch Make/options
echo "LAMMPS_DIR ="$lammps_sir > Make/options

wmake libso dragModels
wmake libso chPressureGrad
wmake libso lammpsFoamTurbulenceModels
wmake
