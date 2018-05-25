#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
#$ -pe orte 5-10
# -pe threaded 5  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
#module load dmtcp
module load R
cd /share/nobackup/noo11/Floc9
R CMD BATCH rnumber.R result1.txt
###############300 design by 10 replicates
#cp AUTO_RESTART_PLEASE k2.sh run/
#cd run
for num in {1..300}
do
layi=$num
echo $layi | sed s/"{1..300}"/$layi/< k2.sh >v$num.sh
#chmod +x v$num.sh
done


for num in {1..300}
do
for fold in {1..10}
do
layi=$fold
echo $layi | sed s/"{1..10}"/$layi/< v$num.sh >v$num.$fold.sh
chmod +x v$num.$fold.sh
done
done
###############submit
for num in {1..2}
do
for fold in {1..3}
do
touch AUTO_RESTART_PLEASE
qsub -P MATH run/v$num.$fold.sh
done
done
###################END

