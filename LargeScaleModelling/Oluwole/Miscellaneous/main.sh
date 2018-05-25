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
module load openmpi-x86_64
module load fftw/3.3.4 
module load R

#cp prep.R  /share/nobackup/noo11/Floc9/inputnum
#cd /share/nobackup/noo11/Floc9/inputnum
#mpirun -np 1 R CMD BATCH prep.R res3.txt
#################Futher preparation
mpirun -np 1 R CMD BATCH prep2.R result2b.txt

