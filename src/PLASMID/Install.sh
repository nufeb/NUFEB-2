# Install/Uninstall package files in LAMMPS
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

# package files without dependencies
action fix_divide_bacillus_minicell.cpp
action fix_divide_bacillus_minicell.h
action fix_property_plasmid.cpp
action fix_property_plasmid.h
action fix_plasmid_conjugation.cpp
action fix_plasmid_conjugation.h
action fix_plasmid_kill.cpp
action fix_plasmid_kill.h
action fix_plasmid_partition.cpp
action fix_plasmid_partition.h
action fix_plasmid_replication.cpp
action fix_plasmid_replication.h
action compute_plasmid_msd.cpp
action compute_plasmid_msd.h
action compute_plasmid_copy.cpp
action compute_plasmid_copy.h
action compute_plasmid_nbirth.cpp
action compute_plasmid_nbirth.h
action compute_plasmid_ave_plasmid.cpp
action compute_plasmid_ave_plasmid.h

if (test $mode = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*PLASMID[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-DLMP_PLASMID |' ../Makefile.package
  fi

elif (test $mode = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*PLASMID[^ \t]* //' ../Makefile.package
  fi

fi