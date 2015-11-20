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
#$ -l s_rt=450:00:00
# Run job through bash shell
#$ -V
#$ -S /bin/bash

cd /frontiers-shared/shared/data/seg-dat/noo11/test3
#########generate random number and LHS with R 
R CMD BATCH rnumber.R result1.txt
##############run the lammp code in parallel


for num in {1..2}
do
cd /frontiers-shared/shared/data/seg-dat/noo11/test3/input$num
#; /frontiers-shared/home/noo11/LAMMPS-Prashant/lammps1Feb2014/src/lmp_serial<Inputscript$num.lammps  
#cd /data/noo11/test3/input$num; /home/noo11/Lammps-NUFEB/lammps1Feb2014/lmp_serial<Inputscript$num.lammps
for fold in {1..5}
do
cd /frontiers-shared/shared/data/seg-dat/noo11/test3/input$num/input$fold; /frontiers-shared/home/noo11/LAMMPS-Prashant/lammps1Feb2014/src/lmp_serial<Inputscript$fold.lammps  
done
cd ../
done
cd ../


###########post-processing with R in parallel
#scp -p -r rlammp.R noo11@segedunum:/frontiers-shared/shared/data/seg-dat/noo11/test2/rlammp2.R
#ssh segedunum
#cd /frontiers-shared/shared/data/seg-dat/noo11/test2
#for num in {1..5}
#do
#done
cd /frontiers-shared/shared/data/seg-dat/noo11/test3; R CMD BATCH rlammp3.R result22.txt  
########################Load data and build the emulator using R
cd /frontiers-shared/shared/data/seg-dat/noo11/test3; R CMD BATCH main2.R result33.txt 



