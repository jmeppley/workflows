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
WORKFLOW_PATH=${3:-/lus/scratch/usr/jmeppley/opt/workflows}

# keep a log of SLURM job IDs in assembly folder
date >> ./.slurm.history
echo -n -e "${SLURM_JOB_ID}\t" >> ./.slurm.history

# report details of job
date
echo "Running job ${SLURM_JOB_ID} in $DIR"
echo -n -e "HOST:\t"
hostname
echo ""
echo loading assembly env...
source activate $CONDA_ENV

COMMAND="snakemake -s${WORKFLOW_PATH}/assembly.metagenomic.snake --configfile ${CONFIG} -p -k -j 20 $DRYRUN --notemp"
echo $COMMAND
echo ""
eval $COMMAND
