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

for fold in 31  45  47  59  74  78  80  81  90  97  98  99 101 103 104 107 108 109 111 126 128 137 142 143 149 153 155 157 167 168 174 178 179 182 183 201 211 213 215 216 222 227 235 237 244 245 250 254 260 266 267 268 276 277 279 285 290 294 298 300
do
cp ana2.R  /share/nobackup/noo11/Floc9/input$fold
cd /share/nobackup/noo11/Floc9/input$fold
ls -d input*/res >  new.txt
R CMD BATCH ana2.R res4.txt
cd ../
done
#################Futher preparation
#mpirun -np 1 R CMD BATCH prep2.R result2b.txt

