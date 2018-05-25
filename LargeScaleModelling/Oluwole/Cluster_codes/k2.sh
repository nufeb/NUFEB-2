#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
# -pe orte 5-20
#$ -pe threaded 4-20  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules

#module load libxml2/2.9.2
export JAVA_HOME=~/jre1.8.0_111
export PATH="$JAVA_HOME/bin:$PATH"
#export JAVA_HOME=~/jdk1.8.0_111
#export PATH="$JAVA_HOME/bin:$PATH"
#cd /share/nobackup/noo11/idyno
#
cd /share/nobackup/noo11/Ca*/idyno2/inputnum
python /home/noo11/iDynoMiCS/sc*/RunIdyno.py /share/nobackup/noo11/Ca*/idyno2/inputnum/multi3D6new*.xml --multiples=5
