#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-405
#SBATCH --output=hiv_tasp.log
#SBATCH --job-name=hiv_tasp
#SBATCH --mem=12G
#SBATCH --time=05:00:00
#SBATCH --mail-type=ALL

# --ntasks-per-node will be used to specify the number
# of cores to use on the machine.

# working directory 
cd /scratch/users/mreitsma/
  
# load modules
ml R

# execute script
Rscript /scratch/users/mreitsma/ppml14/R/simulation_tasp.R $SLURM_ARRAY_TASK_ID "hiv_tasp_low" "run_main.csv" 10 10