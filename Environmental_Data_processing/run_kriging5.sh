#!/bin/bash
#SBATCH -p medium 
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null
#SBATCH --mem=30G
#SBATCH -t 1-00:00:00
module load intel/mkl/64/2017/2.174
echo $PATH
module load openmpi/intel/64/1.10.7
echo $PATH
Rscript 4_run_kriging_all_variables_e5.R 
