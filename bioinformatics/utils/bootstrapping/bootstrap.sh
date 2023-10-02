#!/bin/bash

#SBATCH --job-name=domics # job name
#SBATCH --ntasks=1 # number of nodes to allocate per job
#SBATCH --cpus-per-task=4 # cpus (threads) per task
#SBATCH --partition=standard # partition (queue) to use
#SBATCH --mem 8000 # request 4Gb RAM per node

# example use
# sbatch --array=1-200 ${doomics_repo_path}/bioinformatics/utils/bootstrapping/bootstrap.sh 10000

srun Rscript ${doomics_repo_path}/bioinformatics/utils/bootstrapping/bootstrap.R job_id=${SLURM_ARRAY_TASK_ID} nbs=$1