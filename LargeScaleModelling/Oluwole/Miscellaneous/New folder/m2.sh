#!/bin/bash
#$ -cwd                         # Use current working directory.
# parrallel environment
# -pe orte 10
#$ -pe threaded 20  
# job run time hours:min:sec (always over estimate )
# Run job through bash shell
#$ -S /bin/bash                 # Interpreting script is bash.
#$ -o output.txt                # Recommended (see above comment).
#$ -e error.txt                 # Recommended (see above comment).
# Load required module
module load openmpi-x86_64
module load fftw/3.3.4
module load blas/3.5.0 
module load lapack/r1517
module load R
   
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

#mpirun -np 5 R CMD BATCH test.R result2.txt
R CMD BATCH m22.R m222.txt
fi
exit_status=$?
if [ $exit_status -eq 0 ]
then
#Job has completed; clean up DMTCP files
echo "$(date +%c) ($JOB_ID) : Job completed: cleaning up DMTCP checkpoint files"
rm dmtcp_restart*sh
rm -r ckpt*
fi

