#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-67
#SBATCH --output=scale5.log
#SBATCH --job-name=scale5
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL

# --ntasks-per-node will be used to specify the number
# of cores to use on the machine.

# working directory 
cd /scratch/users/mreitsma/
  
# load modules
ml R

# execute script
Rscript /scratch/users/mreitsma/ppml14/R/simulation.R $SLURM_ARRAY_TASK_ID "scaledrm_0.5" "run_main_scale5.csv" 1 80

