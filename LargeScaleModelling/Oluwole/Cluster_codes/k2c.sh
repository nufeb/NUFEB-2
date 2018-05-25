#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
#$ -pe orte 6-20
# -pe threaded 1-5  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
module load libpng/1.6.19
#module load parallel-studio-xe/2016Update1
module load topsy-openmpi
module load fftw/3.3.4
module load openmpi-x86_64
module load jpeg/9a
module load fftw/2.1.5 
module load openmpi/2.0.1
module load openmpi/2.1.0
module load R
#for num in {11..100}
#do
cp prep1bb.R /share/nobackup/noo11/Ca*/inputnum
cd /share/nobackup/noo11/Ca*/inputnum
R CMD BATCH prep1bb.R res1c.txt
#cd ..
#done
#export PATH="$PATH:$FFTWhome/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$FFTWhome/lib"
#export PATH="$PATH:/usr/lib64/openmpi/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib"

#export PATH="$PATH:$myhome/bin"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$myhome/lib"


#cd /share/nobackup/noo11/review/bio/inputnum/inputfold
#LD_LIBRARY_PATH=/home/noo11/x86_64-linux-gnu /home/noo11/NUFEB/src/lmp_serial -in Inputscript.lammps

~                                                                                  
~                                                   
