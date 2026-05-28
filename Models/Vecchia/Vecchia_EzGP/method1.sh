#!/bin/bash

#SBATCH --account=def-cdlin  # replace this with your supervisors account
#SBATCH --ntasks=1           # number of processes
#SBATCH --mem-per-cpu=8gb     # memory; default unit is megabytes
#SBATCH --time=30:00:00         # time (HH:MM:SS)

module load StdEnv/2020
module load gcc/9.3.0 r/4.2.1



Rscript EzGP_Vecchia_encode_m1.R 
