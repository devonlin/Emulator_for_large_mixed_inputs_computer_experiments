#!/bin/bash

#SBATCH --account=def-cdlin
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=20:00:00
#SBATCH --mem=2gb


module load StdEnv/2020
module load r/4.2.1


#SBATCH --array=1-50
Rscript EzGP_org.R $SLURM_ARRAY_TASK_ID
