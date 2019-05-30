# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files with dependencies

action atom_vec_nufeb.cpp
action atom_vec_nufeb.h
action fix_density.cpp
action fix_density.h
action fix_diffusion_reaction.cpp
action fix_diffusion_reaction.h
action fix_monod.cpp
action fix_monod.h
action fix_monod_het.cpp
action fix_monod_het.h
action grid_vec_monod.cpp
action grid_vec_monod.h
action nufeb_run.cpp
action nufeb_run.h
