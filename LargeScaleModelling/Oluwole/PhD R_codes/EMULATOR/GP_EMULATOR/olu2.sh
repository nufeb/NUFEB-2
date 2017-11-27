
#!/bin/sh
# Your job name
#$ -N olu3
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
export LD_LIBRARY_PATH="/home/physastro/ooyebamiji/lib/:$LD_LIBRARY_PATH"
/usr/bin/Revo64 CMD BATCH had_irri.R EST.txt
