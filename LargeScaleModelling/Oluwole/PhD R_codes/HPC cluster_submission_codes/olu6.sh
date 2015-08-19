#!/bin/sh
# Your job name
#$ -N olu6
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

#cd /padata/alpha/users/ooyebamiji/NEW
/usr/bin/Revo64 CMD BATCH cv_rain_code.R result4.txt


