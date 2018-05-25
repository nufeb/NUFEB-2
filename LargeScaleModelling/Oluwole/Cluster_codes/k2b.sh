#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
#$ -pe orte 1-5
# -pe threaded 1-5
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
#module load libpng/1.6.19
#module load parallel-studio-xe/2016Update1
#module load topsy-openmpi
#module load fftw/3.3.4
#module load openmpi-x86_64
#module load jpeg/9a
#module load fftw/2.1.5 
#module load openmpi/2.0.1
#module load openmpi/2.1.0

#export PATH="$PATH:$FFTWhome/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$FFTWhome/lib"
#export PATH="$PATH:/usr/lib64/openmpi/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib"

#export PATH="$PATH:$myhome/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$myhome/lib"

cd /share/nobackup/noo11/Calibration_inputs/inputnum/inputfold
/home/noo11/NUFEB2.02/nu*/code/lamm*/src/lmp_serial -in Inputscript.lammps
#/home/noo11/NUFEB2.1/nu*/code/lammp*/src/lmp_serial -in Inputscript.lammps
~                                                                                  
~                                                   
