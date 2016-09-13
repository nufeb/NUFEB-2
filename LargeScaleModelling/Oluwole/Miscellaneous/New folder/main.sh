#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
# -pe orte 10
#$ -pe threaded 20
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
#module load mpich-x86_64   
#module load openmpi-x86_64
#module load apps/R
module load openmpi-x86_64
module load fftw/3.3.4
module load blas/3.5.0 
module load lapack/r1517
module load R
   
#LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/frontiers-shared/home/noo11/R/x86_64-redhat-linux-gnu-library/3.2"
#export R_LIBS="/frontiers-shared/home/noo11/R/x86_64-redhat-linux-gnu-library/3.2"
#cd /frontiers-shared/shared/data/seg-dat/noo11/Floc9/floc
#mpirun -np 5 R CMD BATCH test.R result2.txt
 R CMD BATCH new_val2.R result.txt
