#!/bin/sh
# Your job name
#$ -N olu2
# Use current working directory
#$ -cwd
# Join stdout and stderr
#$ -j y
#Run int this project
#$ -P science
# job run time hours:min:sec (always over estimate )
#$ -l s_rt=45:00:00
# Run job through bash shell
#$ -V
#$ -S /bin/bash

#cd /data/noo11/test2 
#########generate random number and LHS with R 
R CMD BATCH rnumber.R result1.txt
##############run the lammp code in parallel
for num in {1..50}
do
cd /data/noo11/test2/input$num; /home/noo11/Lammps-NUFEB/lammps1Feb2014/lmp_serial<Inputscript$num.lammps & 
done

###########post-processing with matlab in parallel
for num in {1..50}
do
cd /data/noo11/train/test2/input$num; matlab -nosplash -nodisplay -nodesktop < /data/noo11/train/test2/readfile.m & 
done
########################Reading excel file into R


