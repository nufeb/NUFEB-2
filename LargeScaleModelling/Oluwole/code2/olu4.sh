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

cd /frontiers-shared/shared/data/seg-dat/noo11/Biofilm2
#########generate random number and LHS with R 
R CMD BATCH rnumber.R result1.txt
##############run the lammp code in parallel
for num in {1..1000}
do
cd /frontiers-shared/shared/data/seg-dat/noo11/Biofilm2/input$num
for fold in {1..100}
do
cd /frontiers-shared/shared/data/seg-dat/noo11/Biofilm2/input$num/input$fold; /frontiers-shared/home/noo11/LAMMPS-Prashant/lammps1Feb2014/src/lmp_serial<Inputscript$fold.lammps &
rm core*
done
cd ../
done
cd ../
#rm -r input*/input*/core*
########################Load data and build the emulator using R
cd /frontiers-shared/shared/data/seg-dat/noo11/Biofilm2; R CMD BATCH main3.R result3.txt 



