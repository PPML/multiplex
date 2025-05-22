#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-101
#SBATCH --output=bind_calib.log
#SBATCH --job-name=bind_calib
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
Rscript /scratch/users/mreitsma/ppml14/R/calibration/bind_runs_calibration.R $SLURM_ARRAY_TASK_ID