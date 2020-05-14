#!/bin/bash
#SBATCH -p medium 
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem=80G
module load intel/mkl/64/2017/2.174
echo $PATH
module load openmpi/intel/64/1.10.7
echo $PATH
Rscript 4_run_kriging_all_variables_e1.R 