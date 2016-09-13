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

cd /share/nobackup/noo11/Floc9/input999
for fold in {1..10}
do
cd input$fold
grep -n "type diameter" *.bubblemd > snap2.txt
wc -l *.bubblemd > tot.txt
cd ../
done


