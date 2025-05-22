#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-18
#SBATCH --output=sens_hcvrr.log
#SBATCH --job-name=sens_hcvrr
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
Rscript /scratch/users/mreitsma/ppml14/R/simulation.R $SLURM_ARRAY_TASK_ID "hcv_rr_tx_on" "run_main_hcv_rr.csv" 1 10


