#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-30
#SBATCH --output=bind_runs.log
#SBATCH --job-name=bind_runs
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL

# --ntasks-per-node will be used to specify the number
# of cores to use on the machine.

# working directory 
cd /scratch/users/mreitsma/
  
# load modules
ml R

# execute script
Rscript /scratch/users/mreitsma/ppml14/R/analysis/bind_runs.R $SLURM_ARRAY_TASK_ID "trans_prob_high"