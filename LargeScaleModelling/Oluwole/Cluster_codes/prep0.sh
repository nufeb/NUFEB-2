#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
#$ -pe orte 1-5
# -pe threaded 3  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
module load openmpi-x86_64
module load fftw/3.3.4
#module load R
#cp prep.R  /share/nobackup/noo11/biofilm/inputnum
#cd /share/nobackup/noo11/biofilm/inputnum
#mpirun -np 1 R CMD BATCH prep.R result2b.txt
#################Futher preparation
#############
for num in {1..10}
do
for fold in {1..5}
do
cd /share/nobackup/noo11/review/bio/input$num/input$fold
grep -n "type diameter" *.bubblemd > snap2
egrep -o '[0-9]+' snap2 > snap2.txt
wc -l *.bubblemd > tot
egrep -o '[0-9]+' tot > tot.txt
#################
#grep "Size" floc > Size
#egrep -o '[0-9]+' Size> Size.txt
#grep "Volume" floc > Volume
#egrep -o '[0-9.e-]+'{6} Volume >  Volume.txt
#grep "Parent floc size" floc > Parent
 #egrep -o '[0-9]+' Parent > Parent.txt
#grep "ntimestep" floc > Timestep
#egrep -o '[0-9]+' Timestep> Timestep.txt
#grep "Average position" floc > Position
# egrep -o '[0-9.e-]+'{7} Position > Position.txt
#grep "floc at ntimestep" floc > Time2
#egrep -o '[0-9]+' Time2> Time2.txt
cd ../
done
#cd ../
done
cd ../
#################

                                                    
