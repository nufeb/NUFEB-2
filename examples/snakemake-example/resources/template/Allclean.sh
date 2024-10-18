#!/usr/bin/env bash 

set -euo pipefail

cd $(dirname $0)

rm snapshot_* || true
rm grid_* || true
rm dump* || true
rm slurm-* || true
rm -rf Results || true
rm output.lammps || true
rm log.lammps || true
rm *.vt* || true
rm done.tkn || true
