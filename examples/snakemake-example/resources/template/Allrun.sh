#!/usr/bin/env bash

set -euo pipefail

cd $(dirname $0)

mpirun -np $2 $1 -in inputscript.nufeb
