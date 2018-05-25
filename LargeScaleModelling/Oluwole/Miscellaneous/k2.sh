#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
#$ -pe orte 2-10
# -pe threaded 3  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required modules
module load dmtcp
module load R
# The following line ensures that the DMTCP coordinator moves with the job when restarted
export DMTCP_HOST=$(hostname)
export DMTCP_PORT=0 
#export TMPDIR=$SGE_CWD_PATH
# Set the checkpoint interval in seconds - adjust accordingly
export DMTCP_CHECKPOINT_INTERVAL=8000
# Launch executable with checkpointing
if [ -a ./dmtcp_restart_script.sh ]
then
# The file dmtcp_restart_script.sh is used to restart the job
echo "$(date +%c) ($JOB_ID) : Restarting job from checkpoint file on $HOSTNAME, job ID $JOB_ID."
./dmtcp_restart_script.sh
else
# This is our first time running the job
echo "$(date +%c) ($JOB_ID) : Running executable with checkpointing for the first time on $HOSTNAME, job ID $JOB_ID."
#dmtcp_launch matlab -nodisplay -nosplash -nojvm < example.m
cd /share/nobackup/noo11/Floc9


#########generate random number and LHS with R
for num in {1..300}
do
cd /share/nobackup/noo11/Floc9/input$num
for fold in {1..10}
do
cd /share/nobackup/noo11/Floc9/input$num/input$fold
dmtcp_launch /home/noo11/new4/NUFEB/src/lmp_serial<Inputscript$fold.lammps
done
cd ../
done
cd ../
##########
fi
exit_status=$?
if [ $exit_status -eq 0 ]
then
#Job has completed; clean up DMTCP files
echo "$(date +%c) ($JOB_ID) : Job completed: cleaning up DMTCP checkpoint files"
rm dmtcp_restart*sh
rm -r ckpt*
fi
