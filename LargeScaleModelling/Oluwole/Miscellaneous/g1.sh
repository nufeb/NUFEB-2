#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
#$ -pe orte 10
# -pe threaded 3  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
module load libpng/1.6.19
module load fftw/3.3.4 
module load openmpi-x86_64
module load jpeg/9a

#export PATH="$PATH:$FFTWhome/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$FFTWhome/lib"
#export PATH="$PATH:/usr/lib64/openmpi/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib"

#export PATH="$PATH:$myhome/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$myhome/lib"

cd /share/nobackup/noo11/Floc9/input72
mpiexec -np 10 /home/noo11/new6/NUFEB/src/lmp_openmpi -in parallel.lammps -p 10x1

