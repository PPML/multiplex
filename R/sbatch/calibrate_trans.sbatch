#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-501
#SBATCH --output=calib_trans.log
#SBATCH --job-name=calib_trans
#SBATCH --mem=12G
#SBATCH --time=01:30:00
#SBATCH --mail-type=ALL

# --ntasks-per-node will be used to specify the number
# of cores to use on the machine.

# working directory 
cd /scratch/users/mreitsma/
  
# load modules
ml R

# execute script
Rscript /scratch/users/mreitsma/ppml14/R/simulation.R $SLURM_ARRAY_TASK_ID "trans_calib" "trans_calibration.csv" 10 10