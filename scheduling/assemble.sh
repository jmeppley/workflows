#!/bin/bash

#SBATCH -n 20
#SBATCH -N 1

#SBATCH -t 14-0
#SBATCH -p scope.q

# %A" is replaced by the job ID
#SBATCH -o assembly.%A.out 
#SBATCH -e assembly.%A.out 

CONFIG=${1:-assembly.yaml}
CONDA_ENV=${2:-assembly}

# keep a log of SLURM job IDs in assembly folder
date >> ./.slurm.history
echo -n -e "${SLURM_JOB_ID}\t" >> ./.slurm.history

# Start log file with details of job
LOG=assembly.log
date >> $LOG
echo "Running job ${SLURM_JOB_ID} in $DIR" >> $LOG
echo -n -e "HOST:\t" >> $LOG
hostname >> $LOG

echo loading assembly env... >> $LOG
source activate $CONDA_ENV

COMMAND="snakemake -s/home/jmeppley/apps/opt/workflows/assembly.metagenomic.snake --configfile ${CONFIG} -p -k -j 20 $DRYRUN --notemp"
echo $COMMAND
echo $COMMAND >> $LOG
echo "" >> $LOG
eval $COMMAND >> $LOG 2>&1
