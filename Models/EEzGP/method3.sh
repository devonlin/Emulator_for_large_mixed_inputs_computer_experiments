#!/bin/bash


#SBATCH --account=def-cdlin  # replace this with your supervisors account
#SBATCH --ntasks=50           # number of processes
#SBATCH --mem-per-cpu=8gb     # memory; default unit is megabytes
#SBATCH --time=10:00:00         # time (HH:MM:SS)

module load StdEnv/2020
module load gcc/9.3.0 r/4.2.1

# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
R -f method3.R

