#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --array=1-150
#SBATCH --output=crn_hiv.log
#SBATCH --job-name=crn_hiv
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL

# --ntasks-per-node will be used to specify the number
# of cores to use on the machine.

# working directory 
cd /scratch/users/mreitsma/
  
# load modules
ml R

# execute script
Rscript /scratch/users/mreitsma/ppml14/R/prep_crn/prep_hiv.R $SLURM_ARRAY_TASK_ID